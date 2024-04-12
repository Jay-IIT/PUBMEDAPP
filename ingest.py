import sys
import logging
import time

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
# connection = None


def get_cursor():
    connection = Connect(**conn_params)
    cursor = connection.cursor()

    return connection, cursor


def update_log_table(processed, error_log=None, status="FINISHED"):
    try:
        global pubmed_log_id
        connection, cursor = get_cursor()
        if error_log:
            cursor.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, ERROR_LOG = %s,
                    STATUS = %s, END_TIME = CURRENT_TIMESTAMP()
                WHERE PUBMED_LOG_ID = %s
            """, (processed, error_log, status, pubmed_log_id))
        else:
            cursor.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, END_TIME = CURRENT_TIMESTAMP(),STATUS = %s
                WHERE PUBMED_LOG_ID = %s
            """, (processed, status, pubmed_log_id))
        cursor.close()
        connection.commit()
    except Exception as e:
        logging.error("Error updating Log Table into Snowflake: %s, for logid %s", e, pubmed_log_id)
        exit(0)


# Function to insert a row into the table
def insert_to_log_table(row):
    try:
        connection, cursor = get_cursor()
        columns = ', '.join(row.keys())
        values = ', '.join('%s' for _ in row)
        sql = f"""
            INSERT INTO BPDWD.AD_HOC.PUBMED_DATA_LOG 
                ({columns},START_TIME)
            VALUES 
                ({values},CURRENT_TIMESTAMP())
        """
        # print(sql)
        # print(tuple(row.values()))
        # Execute the SQL query with values
        cursor.execute(sql, tuple(row.values()))
        cursor.close()
        connection.commit()
    except Exception as e:
        logging.error("Error while inserting log Table into Snowflake: %s, Search key word %s", e,
                      row["SEARCH_KEYWORD"])
        exit(0)


# Function to insert data into the table
def insert_to_datatbl(data_rows):
    try:
        connection, cursor = get_cursor()
        # Extracting data from data_rows and preparing it for bulk insert
        values = [(row['PMID'], row['SEARCH TERM'], row['TITLE'], row['ABSTRACT'],
                   row['AUTHOR'], row['KEYWORDS'], row['PMC'], pubmed_log_id) for row in data_rows]

        # Performing bulk insert
        cursor.executemany("""
            INSERT INTO PUBMED_DATA (PMID, SEARCH_TERM, TITLE, ABSTRACT, AUTHOR_LIST, KEYWORD_LIST, PMCID, PUBMED_LOG_ID)
            VALUES (%s, %s, %s, %s, %s, %s, %s,%s)""", values)
        cursor.close()
        connection.commit()
    except Exception as e:
        logging.error(f"Exception while insert into Data table {e}")
        raise e


def process_batch(args):
    start, end, search_term, search_results, batch_size = args
    logging.info("Going to download record %i to %i" % (start,end))
    try:
        result = pubmed_batch_download(search_term, search_results, batch_size, start)
        if result:
            insert_to_datatbl(result)
            logging.info(f"Downloaded Data {start} to {end},retrieved {len(result)}")
            return len(result)
        logging.error(f"Unable to download {start} to {end}")
        return 0
    except Exception as e:
        logging.error(f"Error inserting data into Snowflake: {e} for Records {start} to {end}")
        raise e


def fetch_and_upload(search_results, search_term, count):
    processed_counts = [0]
    try:
        batch_size = 500
        num_processes = min(10, mp.cpu_count())  # Number of CPU cores
        logging.info(f"Started Multi Processing {time.asctime()}")
        # Define the worker function for processing batches
        def worker(start, end, search_term, search_results, batch_size, processed_counts):
            processed_counts.append(process_batch(start, end, search_term, search_results, batch_size))
        
        # Create a pool of worker processes
        with mp.Pool(processes=num_processes) as pool:
            for start in range(0, count, batch_size):
                end = min(count, start + batch_size)
                pool.apply_async(worker, args=(start, end, search_term, search_results, batch_size, processed_counts))
            
            pool.close()
            pool.join()
        logging.info(f"Ended multi processing {time.asctime()}")
        processed = sum(processed_counts)
        update_log_table(processed)
    except Exception as e:
        logging.error("An error occurred in the fetch_and_upload: %s", e)
        update_log_table(processed, str(e), "Finished with Exceptions")


if __name__ == "__main__":
    try:
        logging.info(f"Ingestion Application Started")
        logging.info(f"pubmed_log_id: {pubmed_log_id}")
        if len(sys.argv) <= 2:
            logging.error("Usage: python ingest.py <search_term>")
            sys.exit(1)
        search_term = sys.argv[1]
        logging.info(f"Received request: {search_term}")
        min_date, max_date = sys.argv[2], sys.argv[3]
        logging.info(f"Received request With Parameters: {search_term},{min_date},{max_date}")
        if min_date == "None" or max_date == "None":
            min_date = None
            max_date = None
        count, search_results = pubmed_search(search_term, min_date=min_date, max_date=max_date)
        log_table_entry = {"PUBMED_LOG_ID": pubmed_log_id,
                           "SEARCH_KEYWORD": search_term,
                           "TOTAL_PMIDS": count,
                           "STATUS": "STARTED",
                           "FROM_DATE": min_date,
                           "TO_DATE": max_date}
        insert_to_log_table(log_table_entry)
        # exit(0)
        fetch_and_upload(search_results, search_term,count)
        logging.info(f"Ingestion Application Stopped Successfully")
    except Exception as e:
        logging.error("An error occurred in the main process: %s", e)
        logging.info(f"Ingestion Application Stopped With Exception")