import streamlit as st
from Bio import Entrez
import re 
import pandas as pd
import numpy as np 

def remove_html_tags(text):
    pattern = re.compile('<.*?>')
    text= re.sub(pattern, '', text)
    text = text.replace("\u2009", " ")
    return text

# Define function to fetch article details from PubMed using Biopython
@st.cache_data 
def fetch_article_details(pubmed_id):
    try:
        # Construct the PubMed query
        handle = Entrez.efetch(db='pubmed', id=pubmed_id)

        # Read and parse the XML response
        record = Entrez.read(handle)
        #print(record)

        # Extract abstracts from the parsed record
        abstracts = []
        try:
            title = record['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']['Title'] 
        except Exception as e:
            title = 'No Title'

        for article in record['PubmedArticle']:
            if 'MedlineCitation' in article:
                citation = article['MedlineCitation']
                if 'Article' in citation:
                    article_info = citation['Article']
                    if 'Abstract' in article_info:
                        abstract = article_info['Abstract']['AbstractText']
                        abstracts.extend(list(map(lambda x : remove_html_tags(x), abstract)))
        if not abstracts:
            abstracts.append("")

        return (pubmed_id,title,abstracts[0])
    except Exception as e:
        st.error(e)


# PubMed search function
@st.cache_data
def search_pubmed(search_term):
    try:
        # Set your email address for Entrez
        Entrez.email = "dummy@yahoo.com"
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=100000)
        record = Entrez.read(handle)
        return record["IdList"]
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
                    data.append({'PMID': result[0],'SEARCH TERM':search_term,'TITLE': result[1], 'ABSTRACT': result[2]})
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
                    search_term = ""
                    st.balloons()
                    st.success("Upload Started [Phase 2]")
        
if __name__ == '__main__':
    try:
        run_streamlit_app()
    except KeyboardInterrupt:
         exit(0)
