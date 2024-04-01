import sys
import logging
from core import search_pubmed_term
from core import fetch_pubmedid_details
from time import sleep
import snowflake.connector

# Set up logging
logging.basicConfig(filename='ingest.log', level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

pkey = None
connection = None
# Establish connection parameters
conn_params = {
    "user": "<your_username>",
    "password": "<your_password>",
    "account": "dn13102.us-east-1",
    "warehouse": "INSIGHTSDWDEVTEST_WH",
    "database": "BPDWD",
    "schema": "AD_HOC"
}

def update_status(pubmed_log_id, processed, error_log=None, status="FINISHED"):
    try:
        cursor = connection.cursor()
        if error_log:
            cursor.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, ERROR_LOG = %s,
                    STATUS = %s, END_TIME = CURRENT_TIMESTAMP()
                WHERE PUBMED_LOG_ID = %s
            """, (processed, error_log, pubmed_log_id, status))
        else:
            cursor.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, END_TIME = CURRENT_TIMESTAMP()
                WHERE PUBMED_LOG_ID = %s
            """, (processed, pubmed_log_id, status))
        cursor.close()
        connection.commit()
    except Exception as e:
        logging.error("Error updating Log Table into Snowflake: %s, for logid %s", e, pubmed_log_id)

# Function to insert a row into the table and return PUBMED_LOG_ID
def insert_row_return_id(search_keyword, total_pmids, status):
    try:
        cursor = connection.cursor()
        cursor.execute("""
            INSERT INTO PUBMED_DATA_LOG 
            (SEARCH_KEYWORD, TOTAL_PMIDS, STATUS, START_TIME)
            VALUES (%s, %s, %s, CURRENT_TIMESTAMP())
            RETURNING PUBMED_LOG_ID
        """, (search_keyword, total_pmids, status))
        pubmed_log_id = cursor.fetchone()[0]
        cursor.close()
        connection.commit()
        return pubmed_log_id
    except Exception as e:
        logging.error("Error while inserting log Table into Snowflake: %s, Search key word %s", e, search_keyword)

# Function to insert data into the table
def insert_to_snowflake(result):
    try:
        # Establish connection
        connection = snowflake.connector.connect(**conn_params)

        # Construct INSERT statement dynamically
        cursor = connection.cursor()
        cursor.execute("""
            INSERT INTO PUBMED_DATA (PMID, SEARCH_TERM, TITLE, ABSTRACT, AUTHOR_LIST,KEYWORD_LIST)
            VALUES (%s, %s, %s, %s, %s, %s)
        """, (result[0], result[1], result[2], result[3], result[4], result[5]))   
        cursor.close()
        connection.commit()

        # Close connection
        connection.close()

    except Exception as e:
        logging.error("Error inserting data into Snowflake: %s", e)

def fetch_and_upload(pubmedids_list):
    processed = 0
    try:
        for idx, pubmedid in enumerate(pubmedids_list, start=1):
            processed += 1
            if idx % 10 == 0:
                sleep(2)  # Introduce a delay every 10 iterations
            result = fetch_pubmedid_details(pubmedid)  # Assuming fetch_pubmedid_details is a function
            if result:
                insert_to_snowflake(result)  # Assuming insert_to_snowflake is a function
            else:
                raise Exception(f"Unable to fetch PubMed ID {pubmedid}")  # Raise exception if fetch failed
        
        update_status(pkey, processed)
    except Exception as e:
        logging.error("An error occurred in the fetch_and_upload: %s", e)
        update_status(pkey, processed, str(e), "Finished with Exceptions")

if __name__ == "__main__":
    try:
        if len(sys.argv) != 2:
            logging.error("Usage: python ingest.py <search_term>")
            sys.exit(1)
        search_term = sys.argv[1]
        pubmedids_list = search_pubmed_term(search_term)
       
        if pubmedids_list:
            connection = snowflake.connector.connect(**conn_params)
            pkey = insert_row_return_id(connection, search_term, len(pubmedids_list), "STARTED")
            fetch_and_upload(pubmedids_list)
    except Exception as e:
        logging.error("An error occurred in the main process: %s", e)
