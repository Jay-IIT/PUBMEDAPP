import streamlit as st

import pandas as pd
import numpy as np 
import subprocess
import os
from core import fetch_pubmedid_details
from core import search_pubmed_term




# Define your search_pubmed, fetch_article_details, and run_streamlit_app functions here...

def run_upload_script(search_term):
    venv_path = os.path.join(os.path.dirname(__file__), 'venv', 'Scripts', 'activate')
    command = f"\"{venv_path}\" && python ingest.py {search_term}"

    subprocess.Popen(command, shell=True, creationflags=subprocess.CREATE_NEW_CONSOLE)


# Define function to fetch article details from PubMed using Biopython
@st.cache_data 
def fetch_article_details(pubmed_id):
    try:
        return fetch_pubmedid_details(pubmed_id)
    except Exception as e:
        st.error(e)
        return None


# PubMed search function
@st.cache_data
def search_pubmed(search_term):
    try:
         return search_pubmed_term(search_term)
    except Exception as e:
        st.error(f"Error occurred during PubMed search: {e}")
        return []
 
    
# Streamlit app code
def run_streamlit_app():

    st.set_page_config(
            page_title="BioPharm Communications Pubmed",
            page_icon="chart_with_upwards_trend",
            layout="wide",
            initial_sidebar_state="collapsed")

    hide_menu_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
    st.markdown(hide_menu_style, unsafe_allow_html=True)
    st.markdown("<h1 style='text-align: center; color: red;'> PUBMED SEARCH </h1>", unsafe_allow_html=True) 
    st.markdown(
            """
        <style>
        button {
            height: auto;
            padding-top: 20px !important;
            padding-bottom: 20px !important;
        }
        </style>
        """,
            unsafe_allow_html=True,
    )
    
    for __ in range(3):
        st.write("\n") 

    # Create three columns layout
    _, col2, col3 = st.columns([1, 3, 1])

    # Add label and text input in the second column
    with col2:
        st.markdown("<div style='font-size: 25px;color: Blue;'>Enter Your Search Term</span></div>", unsafe_allow_html=True)
        search_term = st.text_input(" ", key="search term",value="", help='Type your search term here.')
        search_term = search_term.capitalize()  # Optionally capitalize the input

 

    # Search button to trigger search
    if search_term:
        # Perform search when the search button is clicked
        with st.spinner("Searching..."):
            search_results = search_pubmed(search_term)
            _,col2, col3,col4,_ = st.columns(5)
            reserved = st.container()
            with col2:
                st.markdown(f"<div style='font-size: 20px;'>Total Available PMID'S for {search_term} are  <span style='color: red; font-size: 24px;'><b> {len(search_results)}</b></span></div>", unsafe_allow_html=True)
            if search_results:
                data = []
                num_results_display = 30
                search_results = search_results[:num_results_display]
                with col3:
                    st.markdown(f"<div style='font-size: 20px;'>Displaying Top ↗️📈 <span style='color: red; font-size: 24px;'><b>{len(search_results)}</b></span></div>", unsafe_allow_html=True)
                for pmid in search_results:
                    result = fetch_article_details(pmid)
                    if result:
                        data.append({'PMID': result[0],'SEARCH TERM':search_term,
                                    'TITLE': result[1], 'ABSTRACT': result[2],
                                    'AUTHOR':result[3],'KEYWORDS':result[4]})
                df = pd.DataFrame(data)
                df.index = np.arange(1,len(df)+1)
                st.table(df)
            else:
                st.write("No search results found.")
        if search_results:
            with col4:
                Uploadbtn = st.button("Upload to SnowFlake")
            with reserved:
                if Uploadbtn:
                    st.balloons()
                    run_upload_script(search_term)
                    st.success(f"Upload Started for {search_term}")
                    search_term = ""
        
if __name__ == '__main__':
    try:
        run_streamlit_app()
    except KeyboardInterrupt:
         exit(0)
