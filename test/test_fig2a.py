import pandas as pd

import os

import subprocess

def get_df_fig2a():
    print(os.getcwd())
    fig2a_command = ['./uditas-runner.py', './data/fig2a']

    print(fig2a_command)
    subprocess.call(fig2a_command)

    df = pd.read_excel('./data/fig2a/all_results/fig2a_big_results_pivot.xlsx')

    return df

def get_test_df_fig2a():

    df = pd.read_excel('./test/test_data_all_results/fig2a_big_results_pivot.xlsx')

    return df

def test_fig2a():
    df1 = get_df_fig2a()
    df2 = get_test_df_fig2a()

    assert (df1['Editing Percent'] - df2['Editing Percent']).sum() < 0.00001
