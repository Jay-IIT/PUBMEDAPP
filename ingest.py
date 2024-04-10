import sys
import logging
from core import search_pubmed_term
from core import fetch_pubmedid_details
from time import sleep
import uuid
import multiprocessing as mp
from snowflake.sqlalchemy import URL
from sqlalchemy import create_engine
from snowflake.connector.pandas_tools import pd_writer
import pandas as pd 


# Set up logging
logging.basicConfig(filename='ingest.log', level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

pubmed_log_id = str(uuid.uuid4())
logging.info(f"pubmed_log_id: {pubmed_log_id}")
connection = None
# Establish connection parameters

def create_engine():
    engine = create_engine(URL(user='TEST_USR_ETL',
                               password='8H_3TL@73sT',
                               account='dn13102.us-east-1',
                               warehouse='ETL_TEST',
                               database='BPDWD',
                               schema='AD_HOC'
                               ))
    return engine

def update_log_table(processed, error_log=None, status="FINISHED"):
    try:
        global pubmed_log_id
        engine = create_engine()
        if error_log:
            engine.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, ERROR_LOG = %s,
                    STATUS = %s, END_TIME = CURRENT_TIMESTAMP()
                WHERE PUBMED_LOG_ID = %s
            """, (processed, error_log, status,pubmed_log_id))
        else:
            engine.execute("""
                UPDATE PUBMED_DATA_LOG 
                SET PROCESSED = %s, END_TIME = CURRENT_TIMESTAMP(),STATUS = %s
                WHERE PUBMED_LOG_ID = %s
            """, (processed,status,pubmed_log_id))
        engine.dispose()
    except Exception as e:
        logging.error("Error updating Log Table into Snowflake: %s, for logid %s", e, pubmed_log_id)
        exit(0)

# Function to insert a row into the table  
def insert_to_log_table(row):
    try:
        engine = create_engine()
        df = pd.DataFrame(row)
        df.to_sql("BPDWD.AD_HOC.PUBMED_DATA_LOG", engine, if_exists='replace', index=False,method=pd_writer)
        engine.dispose()
    except Exception as e:
        logging.error("Error while inserting log Table into Snowflake: %s, Search key word %s", e, row["SEARCH_KEYWORD"])
        exit(0)

# Function to insert data into the table
def insert_to_datatbl(data_rows):
    try:
        engine = create_engine()
        df = pd.DataFrame(data_rows)
        df.to_sql("PUBMED_DATA", engine, if_exists='replace', index=False,method=pd_writer)
        engine.dispose()  
    except Exception as e:
        logging.error("Error inserting data into Snowflake: %s", e)
        raise e

def fetch_and_upload(search_term, pubmedids_list):
    processed = 0
    try:
        update_log_table(processed)
    except Exception as e:
        logging.error("An error occurred in the fetch_and_upload: %s", e)
        update_log_table(processed, str(e), "Finished with Exceptions")

if __name__ == "__main__":
    try:
        logging.info(f"Ingestion Application Started")
        if len(sys.argv) != 2:
            logging.error("Usage: python ingest.py <search_term>")
            sys.exit(1)
        search_term = sys.argv[1]
        logging.info(f"Received request: {search_term}")
        pubmedids_list = search_pubmed_term(search_term)
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
                            "FROM_DATE":str(min_date),
                            "TO_DATE":str(min_date),
                            "START_TIME":"CURRENT_TIMESTAMP()"}
        insert_to_log_table(log_table_entry)
        update_log_table(0)
        logging.info(f"Ingestion Application Stopped Successfully")
    except Exception as e:
        logging.error("An error occurred in the main process: %s", e)
        logging.info(f"Ingestion Application Stopped With Exception")
