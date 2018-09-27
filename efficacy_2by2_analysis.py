# Compares the efficacy of two medications through analysis of 2-by-2
# tables

# Copyright (c) 2018 Aubrey Barnard.  This is free software released
# under the MIT License (https://choosealicense.com/licenses/mit/).


# ===== README =====
#
# This program selects data, builds a 2-by-2 table from it, computes
# some statistics on the 2-by-2 table, and then reports all this
# information.  The intention of this program is to provide some
# machinery for comparing two treatments, such as brand and generic
# versions of a medication.
#
# This README uses [Markdown syntax]( https://commonmark.org/) except
# for headings.
#
#
#  ==== Installation ====
#
# There is no need to install this program.  Just place this file
# somewhere that is convenient for your analysis.
#
# The prerequisites are listed in the `import` statements below.  They
# include `barnapy`.  To install `barnapy`, run the following command:
#
#     python3 -m pip install [--user] https://github.com/afbarnard/barnapy/archive/master.zip#egg=barnapy
#
# See https://github.com/afbarnard/fitamord for more information, e.g.,
# whether to use the `--user` option.
#
#
#  ==== Usage ====
#
# To use this program, pass it as the first argument to the Python 3
# interpreter.  Additional command line arguments follow the filename of
# this program, like this:
#
#     python3 <path-to>/efficacy_2by2_analysis.py <sqlite-db> <drug1-id> <drug2-id> <efficacy-sql-file> [<options-file>]
#
# Angle brackets indicate pieces of information that you fill in.
# Square brackets indicate that a piece of information is optional.
#
# See the documentation for `main_cli` for more details on the command
# line interface (CLI) arguments.
#
#
#   === SQLite DB ===
#
# This program assumes you have available a SQLite DB that contains a
# table of drug / prescription records.  Such a table may be created
# from flat files by running [Fitamord](
# https://github.com/afbarnard/fitamord).  Use the options to configure
# this analysis for your DB schema.
#
#
#   === SQL Efficacy Definition ===
#
# The main part of analyzing efficacy is defining what it means for a
# drug to be efficacious and formulating that definition as an SQL query
# (in SQLite dialect: https://sqlite.org/lang.html).  The query must
# return only a single column, and the column must contain subject IDs.
# The query must also include a placeholder ("parameter", denoted as
# `?`) for the drug ID.  This is so the query can be used for both
# drugs.  Once you have formulated your query, save it in a file to pass
# to this program as a CLI argument.
#
# You will probably find it useful to develop and test your query using
# `sqlite3` or another program capable of interfacing with SQLite DBs.
# Note that the efficacy definition query may involve any number of
# tables in the DB, not just the table of drugs.
#
# Here's an example query that defines efficacy as 2 or more drug
# records for a (yet to be determined) drug:
#
#     select person_id
#     from (
#         select person_id, count(*) as n_drug
#         from drug_exposure
#         where drug_concept_id is ?
#         group by person_id
#         having n_drug >= 2
#     );
#
#   === Options ===
#
# The execution of the program can be configured by providing an
# environment file, which contains settings (options) in Python syntax.
# The content between "```" provides an example of what such a file
# might contain.  See the documentation of `default_env` for
# descriptions of individual options.
#
# ```
# # Options for analysis A-12345
# drug_table_name = 'rxs'
# pseudocount = 0.25
# ```
#
#
#  ==== Interpreting the Output ====
#
# The report that this program produces contains two 2-by-2 tables and
# some statistics.  The first table is the expected counts based on the
# observed marginal total counts, and the second table is the observed
# counts.  All tables include the pseudocount in each cell.
#
# Following the tables are some basic statistics of the observed counts.
# Currently, they are the odds ratio (OR) and the relative risk (RR).
# Note that you can copy the counts into your favorite statistical
# program for additional analysis.
#
#
#  ==== Questions ====
#
# Use the source, Luke!  (It should be straightforward and documented.)
# After that, send me e-mail.


import itertools as itools
import math
import random
import sqlite3
import sys

from barnapy import contingency_table as contab
from barnapy import general
from barnapy import parse
#from scipy.stats import poisson


# Environment (dictionary) containing variables that control the
# execution of this analysis.  Edit their values (or provide as a
# command line argument a Python file containing definitions of the same
# variables) to customize this analysis to your situation.
default_env = dict(

    # Name of the table containing prescription records
    drug_table_name='drug_exposure',

    # Name of the column containing the drug ID of the prescription
    drug_id_column_name='drug_concept_id',

    # Name of the column containing the subject / patient ID of the
    # prescription
    subject_id_column_name='person_id',

    # Pseudocount to add to each cell of the 2-by-2 table to smooth it
    pseudocount=0.1,

)


_select_subjects_with_drug_id_sql = 'select {} from {} where {} is ?'

def query_db_to_build_2by2_table(
        db,
        drug1_id,
        drug2_id,
        drug1_efficacy_definition_query_sql,
        drug2_efficacy_definition_query_sql,
        environment,
):
    """
    Query the given DB to build a 2-by-2 table comparing the given
    drugs.

    Return a 2-by-2 table where the exposure is drug 1 (non-exposure is
    drug 2) and the outcome is efficacy according to the given
    definitions.

    db: SQLite DB connection.
    drug1_id, drug2_id: IDs of drug 1 and 2.  Python data type must
        match drug IDs in DB.
    drug{1,2}_efficacy_definition_query_sql: SQL query to identify
        subjects for whom drug 1 (2) is effective.  Must return a
        collection of subject IDs ("select <subject-id> from ...").
        Must have a single parameter, the drug ID ("... where
        ... and <drug-id> is ? ...").
    environment: Dictionary containing definitions describing the table
        of prescriptions: `drug_table_name`, `drug_id_column_name`,
        `subject_id_column_name`.
    """
    # Build query for subjects with a drug
    query = _select_subjects_with_drug_id_sql.format(
        environment['subject_id_column_name'],
        environment['drug_table_name'],
        environment['drug_id_column_name'],
    )
    # Get subjects with drug 1
    cursor = db.execute(query, (drug1_id,))
    subjects1 = set(record[0] for record in cursor)
    # Get subjects with drug 2
    cursor = db.execute(query, (drug2_id,))
    subjects2 = set(record[0] for record in cursor)
    # Get subjects for whom drug 1 is effective
    cursor = db.execute(
        drug1_efficacy_definition_query_sql, (drug1_id,))
    effective1 = set(record[0] for record in cursor)
    # Get subjects for whom drug 2 is effective
    cursor = db.execute(
        drug2_efficacy_definition_query_sql, (drug2_id,))
    effective2 = set(record[0] for record in cursor)
    # Classify the subjects by efficacy.  Trust the sets of subjects
    # rather than the sets of effectives; this allows the efficacy
    # queries to overcount.
    drug1_non = subjects1 - effective1
    drug1_eff = subjects1 - drug1_non
    drug2_non = subjects2 - effective2
    drug2_eff = subjects2 - drug2_non
    # Make the 2-by-2 table.  Exposure is drug1, non-exposure is drug2,
    # outcome is efficacy.
    table = contab.TwoByTwoTable(
        len(drug1_eff), len(drug1_non), # +exp +out | +exp -out
        len(drug2_eff), len(drug2_non), # -exp +out | -exp -out
    )
    return table


def print_table(drug1_id, drug2_id, table, file=sys.stdout):
    header = ('Effective', 'Yes', 'No', 'Totals')
    row_fmt_str = ' {2:<{0}} | {3:<{1}} | {4:<{1}} | {5:<{1}} '
    row_fmt_num = ' {2:>{0}} | {3:>{1}.1f} | {4:>{1}.1f} | {5:>{1}.1f} '
    drug1_id_str = str(drug1_id)
    drug2_id_str = str(drug2_id)
    label_width = max(
        len(header[0]), len(drug1_id_str), len(drug2_id_str))
    table_data = table.table_3x3()
    max_cell = max(itools.chain.from_iterable(table_data))
    cell_width = max(map(len, (*header[1:], '{:.1f}'.format(max_cell))))
    header_str = row_fmt_str.format(label_width, cell_width, *header)
    lines = []
    lines.append(header_str)
    lines.append('-' * len(header_str))
    lines.append(row_fmt_num.format(
        label_width, cell_width, drug1_id_str, *table_data[0]))
    lines.append(row_fmt_num.format(
        label_width, cell_width, drug2_id_str, *table_data[1]))
    lines.append(lines[1])
    lines.append(row_fmt_num.format(
        label_width, cell_width, '', *table_data[2]))
    for line in lines:
        print(line, file=file)


def odds_ratio(table):
    or_ = contab.odds_ratio(table)
    or_se = math.sqrt(1 / table.exp_out + 1 / table.exp_no_out +
                      1 / table.out_no_exp + 1 / table.no_exp_out)
    or_ci_lo = math.exp(math.log(or_) - 1.96 * or_se)
    or_ci_hi = math.exp(math.log(or_) + 1.96 * or_se)
    return or_, or_se, or_ci_lo, or_ci_hi


def relative_risk(table):
    rr = contab.relative_risk(table)
    rr_se = math.sqrt(1 / table.exp_out + 1 / table.exp_tot +
                      1 / table.out_no_exp + 1 / table.no_exp_tot)
    rr_ci_lo = math.exp(math.log(rr) - 1.96 * rr_se)
    rr_ci_hi = math.exp(math.log(rr) + 1.96 * rr_se)
    return rr, rr_se, rr_ci_lo, rr_ci_hi


def expected_counts(table):
    """
    Return a table of expected counts calculated based on the marginal
    counts.
    """
    return contab.TwoByTwoTable(
        table.exp_tot * table.out_tot / table.total,
        table.exp_tot * table.no_out_tot / table.total,
        table.no_exp_tot * table.out_tot / table.total,
        table.no_exp_tot * table.no_out_tot / table.total)


def poisson_test(n_obs, n_exp):
    difference = math.fabs(n_obs - n_exp)
    p_lower = poisson.cdf(n_exp - difference, n_exp)
    p_upper = poisson.sf(n_exp + difference, n_exp)
    return p_lower + p_upper


def poisson_test_table(table):
    # The null hypothesis is expected counts based on marginal counts
    null = expected_counts(table)
    return ((poisson_test(table.exp_out, null.exp_out),
             poisson_test(table.exp_no_out, null.exp_no_out)),
            (poisson_test(table.out_no_exp, null.out_no_exp),
             poisson_test(table.no_exp_out, null.no_exp_out)))


def report_statistics(table, file=sys.stdout):
    or_, _, or_ci_lo, or_ci_hi = odds_ratio(table)
    print('OR: {:<10.3g} \tCI 95%: ({:.3g}, {:.3g})'
          .format(or_, or_ci_lo, or_ci_hi), file=file)
    rr, _, rr_ci_lo, rr_ci_hi = relative_risk(table)
    print('RR: {:<10.3g} \tCI 95%: ({:.3g}, {:.3g})'
          .format(rr, rr_ci_lo, rr_ci_hi), file=file)
    #ptt = poisson_test_table(table)
    #print('Poisson test p-values (2-sided observed vs. expected counts):')
    #print('    {:<9.3g} | {:<9.3g}'.format(*ptt[0]))
    #print('    {:<9.3g} | {:<9.3g}'.format(*ptt[1]))


def main_cli(
        sqlite_db_filename,
        drug1_id_str,
        drug2_id_str,
        efficacy_definition_sql_query_filename,
        environment_filename=None,
        *,
        file=sys.stdout,
):
    """
    Select the counts of the given drugs from the given SQLite DB and
    split them into a 2-by-2 table based on the given definition of
    efficacy.  Report the expected table, observed table, and statistics
    (odds ratio, relative risk) to the given file.

    All the arguments are strings so that this function can handle
    command line arguments.  (Except `file`, which is a text output
    stream like standard output, and is only useful if calling this
    function from another Python program.)

    sqlite_db_filename: Filename of SQLite DB containing table
        of prescription records.
    drug1_id_str: ID of drug 1.
    drug2_id_str: ID of drug 2.
    efficacy_definition_sql_query_filename: Filename of SQL query that
        selects the patients that find a drug efficacious.  The query
        must include a parameter (denoted by `?`) for the drug ID.
    environment_filename: Filename of settings of variables expressed in
        Python syntax.  This is the way to specify additional options
        and configuration.  See the `README` and documentation for
        `default_env`.
    """
    # Parse drug IDs as integers or leave as strings
    drug1_id = parse.int(drug1_id_str, drug1_id_str)
    drug2_id = parse.int(drug2_id_str, drug2_id_str)
    # Assemble efficacy query
    with open(efficacy_definition_sql_query_filename, 'rt') as sql_file:
        efficacy_query = ''.join(sql_file.readlines())
    # Load configuration environment
    environment = dict(default_env)
    if environment_filename is not None:
        environment.update(general.exec_python_file(
            environment_filename))
    # Get counts
    if sqlite_db_filename == 'dummy_data':
        # Generate a dummy 2-by-2 table
        table = contab.TwoByTwoTable(
            *[10 ** random.randrange(5) * random.random()
              for _ in range(4)])
    else:
        # Connect to DB to run queries
        with sqlite3.connect(sqlite_db_filename) as db:
            table = query_db_to_build_2by2_table(
                db,
                drug1_id,
                drug2_id,
                efficacy_query,
                efficacy_query,
                environment,
            )
    obs_counts = table.smoothed(environment['pseudocount'])
    # Report expected counts (smoothed)
    exp_counts = expected_counts(obs_counts)
    print('Expected counts:', file=file)
    print_table(drug1_id, drug2_id, exp_counts, file=file)
    print(file=file)
    # Report observed counts (smoothed)
    print('Observed counts:', file=file)
    print_table(drug1_id, drug2_id, obs_counts, file=file)
    print(file=file)
    # Report statistics
    print('Statistics:', file=file)
    report_statistics(obs_counts, file=file)


if __name__ == '__main__':
    main_cli(*sys.argv[1:])
