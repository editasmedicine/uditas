import pandas as pd

import os

import subprocess

def get_df_fig2b():
    print(os.getcwd())
    fig2b_command = ['./uditas-runner.py', './data/fig2b']

    print(fig2b_command)
    subprocess.call(fig2b_command)

    df = pd.read_excel('./data/fig2b/all_results/fig2b_big_results_pivot.xlsx')

    return df

def get_test_df_fig2b():

    df = pd.read_excel('./test/test_data_all_results/fig2b_big_results_pivot.xlsx')

    return df

def test_fig2b():
    df1 = get_df_fig2b()
    df2 = get_test_df_fig2b()

    assert (df1['Editing Percent'] - df2['Editing Percent']).sum() < 0.00001