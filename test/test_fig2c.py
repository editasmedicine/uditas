import pandas as pd

import os

import subprocess

def get_df_fig2c():
    print(os.getcwd())
    fig2c_command = ['./uditas-runner.py', './data/fig2c']

    print(fig2c_command)
    subprocess.call(fig2c_command)

    df = pd.read_excel('./data/fig2c/all_results/fig2c_big_results_pivot.xlsx')

    return df

def get_test_df_fig2c():

    df = pd.read_excel('./test/test_data_all_results/fig2c_big_results_pivot.xlsx')

    return df

def test_fig2c():
    df1 = get_df_fig2c()
    df2 = get_test_df_fig2c()

    assert (df1['Editing Percent'] - df2['Editing Percent']).sum() < 0.00001