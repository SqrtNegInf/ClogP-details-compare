# ClogP details tables compared side-by-side

The detailed tables for calculated logP are information-dense, and it can be 
difficult to compare results even for members of a congener series.  Placing 
two tables side-by-side can make similarities and differences easier to understand.

## Primary features

* compound input can be SMILES or names
* abbreviations & elisions narrow the combined tables by 20 or more characters
* table rows re-arranged to match 'like' items
* '<' or '>' between tables shows greater-than/less-than relationship
* identical and very similar rows are gray, to 'fade' into background
* annotations to the right of table show direction/magnitude of changes
* 2D images below the table for reference

## Programs

* side-by-side.pl – main routine, input via Web form (CGI), HTML output
* details-filter.pl – helper function, to pre-process raw table

## Caveats

The algorithm will produced a result for any pair of inputs, but the 
output won't be very illuminating unless they share a common scaffold.

## To-do

Producing a more exact concordance between the fragments present in the pair
of compounds would allow a better alignment of Proximity corrections, but it 
would be significantly more work (and never reach 100% matching).
