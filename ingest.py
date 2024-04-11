import sys
import logging
from core import pubmed_search
from core import pubmed_batch_download
import uuid
import multiprocessing as mp
from snowflake.connector import connect as Connect
 


# Set up logging
logging.basicConfig(filename='ingest.log', level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

conn_params = {
    "user": "TEST_USR_ETL",
    "password": "8H_3TL@73sT",
    "account": "dn13102.us-east-1",
    "warehouse": "ETL_TEST",
    "database": "BPDWD",
    "schema": "AD_HOC"
}


pubmed_log_id = str(uuid.uuid4())
logging.info(f"pubmed_log_id: {pubmed_log_id}")
connection = None
 

def get_cursor():
    connection = Connect(**conn_params)
    cursor = connection.cursor()

    return connection,cursor

def update_log_table(processed, error_log=None, status="FINISHED"):
    try:
        global pubmed_log_id
        cursor,connection = get_cursor()
        if error_log:
            cursor.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, ERROR_LOG = %s,
                    STATUS = %s, END_TIME = CURRENT_TIMESTAMP()
                WHERE PUBMED_LOG_ID = %s
            """, (processed, error_log, status,pubmed_log_id))
        else:
            cursor.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, END_TIME = CURRENT_TIMESTAMP(),STATUS = %s
                WHERE PUBMED_LOG_ID = %s
            """, (processed,status,pubmed_log_id))
        cursor.close()
        connection.commit()
    except Exception as e:
        logging.error("Error updating Log Table into Snowflake: %s, for logid %s", e, pubmed_log_id)
        exit(0)

# Function to insert a row into the table  
def insert_to_log_table(row):
    try:
        cursor,connection = get_cursor()
        columns = ', '.join(row.keys())
        values = ', '.join('%s' for _ in row)
        sql = f"""
            INSERT INTO BPDWD.AD_HOC.PUBMED_DATA_LOG 
                ({columns})
            VALUES 
                ({values})
        """

        # Execute the SQL query with values
        cursor.execute(sql, tuple(row.values()))
        cursor.close()
        connection.commit() 
    except Exception as e:
        logging.error("Error while inserting log Table into Snowflake: %s, Search key word %s", e, row["SEARCH_KEYWORD"])
        exit(0)

# Function to insert data into the table
def insert_to_datatbl(data_rows, cursor):
    try:
        cursor,connection = get_cursor()
        # Extracting data from data_rows and preparing it for bulk insert
        values = [(row['PMID'], row['SEARCH_TERM'], row['TITLE'], row['ABSTRACT'],
                   row['AUTHOR_LIST'], row['KEYWORD_LIST'], row['PMCID'] ,pubmed_log_id) for row in data_rows]

        # Performing bulk insert
        cursor.executemany("""
            INSERT INTO PUBMED_DATA (PMID, SEARCH_TERM, TITLE, ABSTRACT, AUTHOR_LIST, KEYWORD_LIST,PMCID,PUBMED_LOG_ID)
            VALUES (%s, %s, %s, %s, %s, %s, %s,%s)""", values)
        cursor.close()
        connection.commit() 
    except Exception as e:
        raise e

def process_batch(args):
    start, end, search_term, search_results, batch_size = args
    try:
        result = pubmed_batch_download(search_term, search_results, batch_size, start)
        insert_to_datatbl(result)
        logging.info(f"Downloaded Data {start} to {end}")
        return end - start  
    except Exception as e:
        logging.error("Error inserting data into Snowflake: {e} for Records {start} to {end}")
        raise e
         

def fetch_and_upload(search_results, search_term):
    try:
        batch_size = 500
        num_processes = min(10,mp.cpu_count()) # Number of CPU cores

        # Create a process pool
        with mp.Pool(processes=num_processes) as pool:
            args_list = [(start, min(count, start + batch_size), search_term, search_results, batch_size) 
                            for start in range(0, count, batch_size)]
            processed_counts = pool.map(process_batch, args_list)

        # After all processes are finished, update log table
        processed = sum(processed_counts) + 1
        update_log_table(processed)
    except Exception as e:
        logging.error("An error occurred in the fetch_and_upload: %s", e)
        update_log_table(processed, str(e), "Finished with Exceptions")

if __name__ == "__main__":
    try:
        logging.info(f"Ingestion Application Started")
        if len(sys.argv) <= 2:
            logging.error("Usage: python ingest.py <search_term>")
            sys.exit(1)
        search_term = sys.argv[1]
        logging.info(f"Received request: {search_term}")
        min_date,max_date = sys.argv[2],sys.argv[3]
        logging.info(f"Received request With Parameters: {search_term},{min_date},{max_date}")
        if min_date == "None" or max_date == "None":
            min_date = None
            max_date = None
        count,search_results = pubmed_search(search_term,min_date=None,max_date=None)
        log_table_entry = {"PUBMED_LOG_ID":pubmed_log_id,
                           "SEARCH_KEYWORD":search_term, 
                            "TOTAL_PMIDS":count, 
                            "STATUS":"STARTED",
                            "FROM_DATE":str(min_date) if min_date is not None else "N/A",
                            "TO_DATE":str(max_date) if max_date is not None else "N/A",
                            "START_TIME":"CURRENT_TIMESTAMP()"}
        insert_to_log_table(log_table_entry)
        fetch_and_upload(search_results,search_term)
        logging.info(f"Ingestion Application Stopped Successfully")
    except Exception as e:
        logging.error("An error occurred in the main process: %s", e)
        logging.info(f"Ingestion Application Stopped With Exception")
