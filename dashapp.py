import pandas as pd
import dash
from dash.dependencies import Input, Output
from dash import dash_table
from dash import html
 
 
 
def run_dash_app(df):
    global dash_app
    dash_app = dash.Dash(__name__)    
    PAGE_SIZE = 5    
    
    dash_app.layout = html.Div([    
        html.Div(dash_table.DataTable(    
        id='datatable-paging',    
        columns=[{"name": i, "id": i} for i in df.columns],  # Exclude index column
        data=df.to_dict(orient='records'),
        page_size=PAGE_SIZE,
        page_current=0,
        filter_action='native',   # Enable custom filtering
        sort_action='native', 
        style_cell_conditional=[
            {'if': {'column_id': c},
            'whiteSpace': 'pre-line',  # Set whiteSpace to pre-line for multiline text
            } for c in df.columns
        ],
        style_filter_conditional=[
                    {'if': {'column_id': c},
                    'textAlign': 'left',  # Left-align filter data
                    } for c in df.columns
        ],
        style_cell={'textAlign': 'left'},  # Left-align text in cells
        style_filter={'textAlign': 'left'},  # Left-align filter text
        style_header={'textAlign': 'left'},  # Left-align header
        style_table={'marginLeft': '10px'}  # Adjust left margin  # Adjust margins
    ))    
    ])        

    @dash_app.callback(    
        Output('datatable-paging', 'page_size'),    
        [Input('select_page_size', 'value')])    
    def update_graph(page_size):    
        return page_size
    


    try:
        dash_app.run_server(debug=False)
    except KeyboardInterrupt:
        dash_app.server.stop()
    