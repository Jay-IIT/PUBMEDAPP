import streamlit as st
from dashapp import run_dash_app
import threading
from Bio import Entrez
import re 
import pandas as pd

dash_app = None

def stop_dash_app():
    global dash_app
    if dash_app:
        dash_app.server.stop()



def remove_html_tags(text):
    pattern = re.compile('<.*?>')
    text= re.sub(pattern, '', text)
    text = text.replace("\u2009", " ")
    return text

# Define function to fetch article details from PubMed using Biopython
def fetch_article_details(pubmed_id):
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

    return (pubmed_id,title,abstracts)


# PubMed search function
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
            padding-top: 10px !important;
            padding-bottom: 10px !important;
        }
        </style>
        """,
            unsafe_allow_html=True,
    )
    
    for __ in range(3):
        st.write("\n") 


     # Autocomplete text input for entering search terms
    search_term = st.text_input("Enter Pubmed Search term")
    search_term = search_term.capitalize()

    

    # Search button to trigger search
    if st.button("Search",):
        # Perform search when the search button is clicked
        with st.spinner("Searching..."):
            search_results = search_pubmed(search_term)
            st.markdown(f"<div style='font-size: 20px;'>Total PMID'S available for Search Term {search_term} are  <span style='color: red; font-size: 24px;'><b>{len(search_results)}</b></span></div>", unsafe_allow_html=True)
            if search_results:
                data = []
                for pmid in search_results[:30]:
                    result = fetch_article_details(pmid)
                    data.append({'ID': result[0], 'TITLE': result[1], 'ABSTRACT': result[2]})
                df = pd.DataFrame(data)
                stop_dash_app()
                dash_thread = threading.Thread(target=run_dash_app,args=(df,),daemon=True)
                dash_thread.start() 
                st.components.v1.iframe("http://localhost:8050", width=1500, height=2000)
               # st.markdown(
               #       """
               #      <div style="display: flex; justify-content: center;">
               #          <iframe src="http://localhost:8050" width="1000" height="1000"></iframe>
               #      </div>
               #      """,
               #      unsafe_allow_html=True
               #  )
            else:
                st.write("No search results found.")
        


   #st.components.v1.iframe("http://localhost:8050", width=1000, height=1000)
   # st.markdown(
   #     """
   #    <div style="display: flex; justify-content: center;">
   #        <iframe src="http://localhost:8050" width="1000" height="1000"></iframe>
   #    </div>
   #    """,
   #    unsafe_allow_html=True
   #)


if __name__ == '__main__':
   # dash_thread = threading.Thread(target=run_dash_app,daemon=True)
   # dash_thread.start()
    try:
        run_streamlit_app()
    except KeyboardInterrupt:
        stop_dash_app()
