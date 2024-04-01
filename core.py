from Bio import Entrez
import re 


def fetch_pubmedid_details(pubmed_id):
    def remove_html_tags(text):
        try:
            pattern = re.compile('<.*?>')
            text= re.sub(pattern, '', text)
            text = text.replace("\u2009", " ")
            return text
        except Exception as e:
            return " "
    try:
        # Construct the PubMed query
        handle = Entrez.efetch(db='pubmed', id=pubmed_id)

        # Read and parse the XML response
        record = Entrez.read(handle)

        # Extract abstracts from the parsed record
        abstracts = []
        try:
            title = record['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']['Title'] 
        except Exception as e:
            title = 'No Title'
        try :
            authors = [author["ForeName"] +" "+author["LastName"]
                       for author in record['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']]
        except Exception as e:
            authors = []   
        try:
            keywords = [kw for kw in record['PubmedArticle'][0]['MedlineCitation']["KeywordList"][0]]
        except Exception as e:
            keywords = []
        authors = " ,".join(authors)
        keywords = " ,".join(keywords)
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

        return (pubmed_id,title,abstracts[0],authors,keywords)
    except Exception as e:
        return None

def search_pubmed_term(search_term):
    try:
        # Set your email address for Entrez
        Entrez.email = "dummy@yahoo.com"
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=100000)
        record = Entrez.read(handle)
        return record["IdList"]
    except Exception as e:
        return []
    
