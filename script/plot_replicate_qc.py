import plotly.graph_objs as go
import plotly.offline as pyo
import pandas as pd


df = pd.read_csv('result/mut_scores.tsv', sep='\t')

correlation = df['Rep1_exp_score'].corr(df['Rep2_exp_score'])

hovertext = ['<br>'.join(['SH_pep: ' + str(shpep),
                  'S2P6 mut: ' + str(s2p6_mut),
                  'Rep1 expression: ' + str(rep1_enrich),
                  'Rep2 expression: ' + str(rep2_enrich),]
                    )
      for shpep, s2p6_mut, rep1_enrich, rep2_enrich in zip(df['SH_pep'], df['mut_ID'], df['Rep1_exp_score'],df['Rep2_exp_score'])]
trace = go.Scatter(x=df['Rep1_exp_score'], y=df['Rep2_exp_score'], mode='markers', hovertext=hovertext)
data = [trace]
layout = go.Layout(title='Expression correlation',
            xaxis=dict(title='Rep1 expression'),
            yaxis=dict(title='Rep2 expression'),
            annotations=[
        # Add the correlation coefficient as a text annotation
          dict(
            text=f'Correlation: {correlation:.2f}',
            xref='paper', yref='paper',  # Use relative coordinates
            x=0.95, y=0.05,  # Adjust these values for positioning
            showarrow=False,
            font=dict(size=12, color='black'),
            bgcolor='lightgray',
            bordercolor='gray',
            borderwidth=2,
            borderpad=4,
            opacity=0.8,
        )
    ])
fig = go.Figure(data=data, layout=layout)

pyo.plot(fig, filename='graph/QC/Expression_score_correlation.html')


correlation = df['Rep1_binding_score'].corr(df['Rep2_binding_score'])

hovertext = ['<br>'.join(['SH_pep: ' + str(shpep),
                  'S2P6 mut: ' + str(s2p6_mut),
                  'Rep1 expression: ' + str(rep1_enrich),
                  'Rep2 expression: ' + str(rep2_enrich),]
                    )
      for shpep, s2p6_mut, rep1_enrich, rep2_enrich in zip(df['SH_pep'], df['mut_ID'], df['Rep1_binding_score'],df['Rep2_binding_score'])]
trace = go.Scatter(x=df['Rep1_binding_score'], y=df['Rep2_binding_score'], mode='markers', hovertext=hovertext)
data = [trace]
layout = go.Layout(title='Binding correlation',
            xaxis=dict(title='Rep1 binding'),
            yaxis=dict(title='Rep2 binding'),
            annotations=[
        # Add the correlation coefficient as a text annotation
          dict(
            text=f'Correlation: {correlation:.2f}',
            xref='paper', yref='paper',  # Use relative coordinates
            x=0.95, y=0.05,  # Adjust these values for positioning
            showarrow=False,
            font=dict(size=12, color='black'),
            bgcolor='lightgray',
            bordercolor='gray',
            borderwidth=2,
            borderpad=4,
            opacity=0.8,
        )
    ])
fig = go.Figure(data=data, layout=layout)

pyo.plot(fig, filename='graph/QC/Binding_score_correlation_1*10-5.html')
