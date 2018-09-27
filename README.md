Analysis of Drug Efficacy via 2-by-2 Tables
===========================================

This is a one-off repository to deliver a Python script to the world.
As such, I took the repository name GitHub suggested.  :grinning:
:octocat:

The program, `efficacy_2by2_analysis.py`, is a Python 3 program that
selects data from a SQLite DB, builds a 2-by-2 table from it, computes
some statistics on the 2-by-2 table, and then reports all this
information.  The intention of this program is to provide some machinery
for comparing two treatments, such as brand and generic versions of a
medication.


Installation
------------

There is no need to install this program.  Just place the script
somewhere that is convenient for your analysis.

The only notable prerequisite is `barnapy`.  To install `barnapy`, run
the following command:

    python3 -m pip install [--user] https://github.com/afbarnard/barnapy/archive/master.zip#egg=barnapy

See https://github.com/afbarnard/fitamord for more information, e.g.,
whether to use the `--user` option.


Usage
-----

To use this program, pass it as the first argument to the Python 3
interpreter.  Additional command line arguments follow the filename of
this program, like this:

    python3 <path-to>/efficacy_2by2_analysis.py <sqlite-db> <drug1-id> <drug2-id> <efficacy-sql-file> [<options-file>]

Angle brackets indicate pieces of information that you fill in.  Square
brackets indicate that a piece of information is optional.

See the documentation for `main_cli` for more details on the command
line interface (CLI) arguments.


### SQLite DB ###

This program assumes you have available a SQLite DB that contains a
table of drug / prescription records.  Such a table may be created from
flat files by running [Fitamord](
https://github.com/afbarnard/fitamord).  Use the options to configure
this analysis for your DB schema.


### SQL Efficacy Definition ###

The main part of analyzing efficacy is defining what it means for a drug
to be efficacious and formulating that definition as an SQL query (in
SQLite dialect: https://sqlite.org/lang.html).  The query must return
only a single column, and the column must contain subject IDs.  The
query must also include a placeholder ("parameter", denoted as `?`) for
the drug ID.  This is so the query can be used for both drugs.  Once you
have formulated your query, save it in a file to pass to this program as
a CLI argument.

You will probably find it useful to develop and test your query using
`sqlite3` or another program capable of interfacing with SQLite DBs.
Note that the efficacy definition query may involve any number of tables
in the DB, not just the table of drugs.

Here's an example query that defines efficacy as 2 or more drug records
for a (yet to be determined) drug:

    select person_id
    from (
        select person_id, count(*) as n_drug
        from drug_exposure
        where drug_concept_id is ?
        group by person_id
        having n_drug >= 2
    );


### Options ###

The execution of the program can be configured by providing an
environment file, which contains settings (options) in Python syntax.
The following content provides an example of what such a file might
contain.  See the documentation of `default_env` for descriptions of
individual options.

```
# Options for analysis A-12345
drug_table_name = 'rxs'
pseudocount = 0.25
```


Interpreting the Output
-----------------------

The report that this program produces contains two 2-by-2 tables and
some statistics.  The first table is the expected counts under the null
hypothesis that the exposure and outcome are independent based on the
observed marginal total counts.  The second table is the observed
counts.  All tables include the pseudocount in each cell.

Following the tables are some basic statistics of the observed counts.
Currently, they are the odds ratio (OR) and the relative risk (RR).
Note that you can copy the counts into your favorite statistical
program for additional analysis.


-----

Copyright (c) 2018 Aubrey Barnard.  This is free software released under
the MIT License.  See `LICENSE` for details.
