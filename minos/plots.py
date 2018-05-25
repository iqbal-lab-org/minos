import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
#from matplotlib.lines import Line2D


def load_dp_and_gt_conf_data_from_file(infile):
    return pd.read_csv(infile, header=0, sep='\t')


def scatter_plot_gt_conf_vs_dp(data, outfile):
    p = sns.lmplot(
        x="DP",
        y="GT_CONF",
        data=data,
        fit_reg=False,
        legend=True,
        markers=["o"],
        scatter_kws={"s": 5},
    )
    p.savefig(outfile)
    plt.clf()


def scatter_plot_gt_conf_vs_dp_colour_by_tp_fp(data, outfile):
    p = sns.lmplot(
        x="DP",
        y="GT_CONF",
        data=data,
        fit_reg=False,
        hue='TP_OR_FP',
        legend=True,
        markers=["o", "x"],
        scatter_kws={"s": 5},
        palette=dict(TP="g", FP="m")
    )
    p.savefig(outfile)
    plt.clf()


def histogram_of_one_dataframe_column(data, column_name, outfile):
    p = sns.distplot(data[column_name], kde=False)
    p.figure.savefig(outfile)
    plt.clf()


def histogram_of_one_dataframe_column_color_by_tp_fp(data, column_name, outfile):
    tps = data[data['TP_OR_FP'] == 'TP']
    fps = data[data['TP_OR_FP'] == 'FP']
    p = sns.distplot(tps[column_name], color="g", kde=False)
    p = sns.distplot(fps[column_name], color="m", kde=False)
    p.figure.savefig(outfile)
    plt.clf()

