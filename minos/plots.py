import logging

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')

from cluster_vcf_records import vcf_file_read

def load_dp_and_gt_conf_data_from_file(infile):
    return pd.read_csv(infile, header=0, sep='\t')


def scatter_plot_gt_conf_vs_dp(data, outfile):
    p = sns.lmplot(
        x="DP",
        y="GT_CONF",
        data=data,
        fit_reg=False,
        legend=False,
        markers=["o"],
        scatter_kws={"s": 5},
    )
    p.savefig(outfile)
    plt.clf()


def scatter_plot_gt_conf_vs_dp_colour_by_tp_fp(data, outfile, tp_or_fp_types):
    colours = {'TP': 'g', 'FP':'m', 'Unknown': 'black'}
    palette = {x: colours[x] for x in tp_or_fp_types}

    p = sns.lmplot(
        x="DP",
        y="GT_CONF",
        data=data,
        fit_reg=False,
        hue='TP_OR_FP',
        legend=True,
        legend_out=True,
        markers=["o"] * len(palette),
        scatter_kws={"s": 5},
        palette=palette
    )
    p._legend.set_title(None)
    p.savefig(outfile)
    plt.clf()


def histogram_of_one_dataframe_column(data, column_name, outfile):
    p = sns.distplot(data[column_name], kde=False)
    p.figure.savefig(outfile)
    plt.clf()


def histogram_of_one_dataframe_column_color_by_tp_fp(data, column_name, outfile):
    tps = data[data['TP_OR_FP'] == 'TP']
    fps = data[data['TP_OR_FP'] == 'FP']
    p = sns.distplot(tps[column_name], label='True positive', color="g", kde=False)
    p = sns.distplot(fps[column_name], label='False positive', color="m", kde=False)
    plt.legend()
    p.figure.savefig(outfile)
    plt.clf()


def minos_vcf_to_plot_data(infile, outfile):
    header, records = vcf_file_read.vcf_file_to_list(infile)
    data = []
    tp_or_fp_types = set()
    output_cols = ['DP', 'GT_CONF']

    for record in records:
        if record.FORMAT is None:
            continue

        check_geno = record.FORMAT.get('MINOS_CHECK_GENOTYPE', None)
        dp = record.FORMAT.get('DP', None)
        gt_conf = record.FORMAT.get('GT_CONF', None)
        if dp is not None and gt_conf is not None:
            to_append = [dp, gt_conf]

            if check_geno is not None:
                if check_geno == '0':
                    tp_or_fp_type = 'FP'
                elif check_geno == '1':
                    tp_or_fp_type = 'TP'
                else:
                    tp_or_fp_type = 'Unknown'

                tp_or_fp_types.add(tp_or_fp_type)
                to_append.append(tp_or_fp_type)

            data.append(to_append)

    if len(data) == 0:
        logging.warning('No DP and GT_CONF data found in VCF file ' + infile + ' therefore no plots will be made')
        return None

    if 'TP' in tp_or_fp_types or 'FP' in tp_or_fp_types:
        output_cols.append('TP_OR_FP')

    with open(outfile, 'w') as f:
        print(*output_cols, sep='\t', file=f)
        for l in data:
            if len(l) < len(output_cols):
                l.append('Unknown')
            print(*l[:len(output_cols)], sep='\t', file=f)

    return tp_or_fp_types


def plots_from_minos_vcf(infile, outprefix):
    data_tsv = outprefix + '.data.tsv'
    tp_or_fp_types = minos_vcf_to_plot_data(infile, data_tsv)
    if tp_or_fp_types is None:
        return

    data = load_dp_and_gt_conf_data_from_file(data_tsv)
    scatter_file = outprefix + '.gt_conf_dp_scatter.pdf'
    dp_hist_file = outprefix + '.dp_hist.pdf'
    gt_conf_hist_file = outprefix + '.gt_conf_hist.pdf'

    if len(tp_or_fp_types) == 0 or tp_or_fp_types == {'Unknown'}:
        scatter_plot_gt_conf_vs_dp(data, scatter_file)
        histogram_of_one_dataframe_column(data, 'DP', dp_hist_file)
        histogram_of_one_dataframe_column(data, 'GT_CONF', gt_conf_hist_file)
    else:
        scatter_plot_gt_conf_vs_dp_colour_by_tp_fp(data, scatter_file, tp_or_fp_types)
        histogram_of_one_dataframe_column_color_by_tp_fp(data, 'DP', dp_hist_file)
        histogram_of_one_dataframe_column_color_by_tp_fp(data, 'GT_CONF', gt_conf_hist_file)

