from Bio import Entrez
import re 
import xmltodict


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
    

def pubmed_search(search_term,min_date=None,max_date=None):
    try:
        # Set your email address for Entrez
        Entrez.email = "dummy@yahoo.com"
        search_results = Entrez.read(
            Entrez.esearch(
                db="pubmed", term=search_term, datetype="pdat", usehistory="y",
                mindate=min_date,maxdate=max_date
            )
        )
        count = int(search_results["Count"]) 
        return count,search_results
    except Exception as e:
        return []

def get_details(article,search_term):
    parsed_result = dict()
    parsed_result['SEARCH TERM'] = search_term
    try:
        parsed_result['TITLE'] = article['MedlineCitation']['Article']['Journal']['Title']
    except:
        parsed_result['TITLE'] = " "
    try:
        parsed_result['PMID'] = article['MedlineCitation']['PMID']["#text"]
    except:
        parsed_result['PMID'] = " "

    try:
        parsed_result['ABSTRACT'] = article['MedlineCitation']['Article']['Abstract']['AbstractText']["#text"]
    except:
        parsed_result['ABSTRACT'] = None
    
    if parsed_result['ABSTRACT'] is None:
        try:
            parsed_result['ABSTRACT'] = " "
            for abstract in article['MedlineCitation']['Article']['Abstract']['AbstractText']:
                parsed_result['ABSTRACT'] += abstract["#text"]
        except Exception as e:
                parsed_result['ABSTRACT'] = None

    if parsed_result['ABSTRACT'] is None:
        try:
            parsed_result['ABSTRACT'] =  article['MedlineCitation']['Article']['Abstract'] 
        except Exception as e:
                parsed_result['ABSTRACT'] = "  "

    if "AbstractText" in parsed_result['ABSTRACT'] :
        try:
            parsed_result['ABSTRACT'] =  parsed_result['ABSTRACT']["AbstractText"]
        except Exception as e:
                parsed_result['ABSTRACT'] = "  "         

    try:
        parsed_result['AUTHOR'] = ",".join([author["LastName"]+" "+author["ForeName"] 
                                            for author in article['MedlineCitation']['Article']['AuthorList']['Author']])
    except:
        parsed_result['AUTHOR'] = " "

    try:
        parsed_result['KEYWORDS'] = ",".join([keyword["#text"] 
                                            for keyword in article['MedlineCitation']["KeywordList"]["Keyword"]])
    except:
        parsed_result['KEYWORDS'] = " "
    
    return parsed_result
    
def pubmed_batch_download(search_term,search_results,batch_size,start=0):
    try:
        stream = Entrez.efetch(
            db="pubmed",
            rettype="medline",
            retmode="xml",
            retstart=start,
            retmax=batch_size,
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
        )
        data = stream.read()
        stream.close()
        # Parse the XML data into a dictionary
        parsed_data = xmltodict.parse(data)
        result = []
        try:
            articles = parsed_data["PubmedArticleSet"]['PubmedArticle']
        except:
            return result

        if type(articles) == dict():
            return get_details(articles,search_term)
        for article in articles:
            result.append(get_details(article,search_term))
         
            
        return result
    except Exception as e:
        return []
    