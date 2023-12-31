﻿# WDS_Gaia

Important information for navigating Daphne's github repository

The 'pipeline' for this project is spread among a few different files. This was my first programming project so it's not the most efficient code, but I tried to leave comments to explain what I was doing.

Here's the 'trajectory' for the project.

1) Prepare the WDS Table: 06-16-22_build_and_write_astropy_table.py
    - this code reads the WDS text file, turns it into a table, and prepares the table to be in a useful format for querying Gaia.
        - generates the secondary coordinate from the primary coordinate, offset, and separation angle
        - updates the coordinates using the proper motion (J200 -> J2016 if the object has proper motion information)

2) Query Gaia: final_query_3-23.ipynb
    - this code divides the WDS table and queries Gaia
    - it takes a long time to run this, so the more cores (to parallelize) the better

3) Fix The Bugs in the Query Results Table: fix_query_results_table4-19.ipynb
    - There were some duplicate rows that needed to be removed
    - Another book-keeping error: the code didn't properly handle some of the entries where Gaia couldn't find 2 stars for a system. These systems were removed from the query results table and added to the table of what I call "index error queries."

4) Sort Query Results Based on Parallax, Separation, and Proper Motion: sort_query_results_4_11.ipynb
    - Use the methods from El-Badry to sort the systems
  

Finally, publication_plots_9-2-23.ipynb is the notebook where I've made all of the plots
