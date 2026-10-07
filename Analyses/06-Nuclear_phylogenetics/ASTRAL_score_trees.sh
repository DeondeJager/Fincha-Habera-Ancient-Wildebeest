# Score astral extended species tree after initial run (if no final normalized quartet score is provided "because of the existense of polytomies and to save time").

## The above message will be in the log file when a mapping (population/species) file is used, and at the end of the log file will be the "extended species tree" (containing individuals) you need to score to get the normalized quartet score - copy and paste it into a new file: e.g. n82_fossil.50kb.ASTRAL.spTree

## But this tree also has additional info (i.e. posterior probability for all three alternatives) which makes it so that astral itself can't read it due to its notation '[pp1=x;pp2=y;pp3=z]' and throws an error ("expected a ')'"). So you have to manually remove these annotations in a text editor and then score the tree:

## ASTRAL annotation options with -t:
#  [(-t|--branch-annotate) <branch annotation level>]
#        How much annotations should be added to each branch: 0, 1, or 2.
#        0: no annotations.
#        1: only the quartet support for the main resolution.
#        2: full annotation (quartet support, quartet frequency, and posterior
#        probability for all three alternatives, plus total number of quartets
#        around the branch and effective number of genes).
#        3 (default): only the posterior probability for the main resolution.
#        4: three alternative posterior probabilities.
#        8: three alternative quartet scores.
#        16/32: hidden commands useful to create a file called freqQuad.csv.
#        10: p-values of a polytomy null hypothesis test. (default: 3)

# Score tree with -t 1 (as instructed by the log file) to get the final quartet score for the overall tree and for each branch:
ASTRAL=/home/pzx702/programs/Astral
java -jar $ASTRAL/astral.5.7.8.jar -i n82_fossil.wOutgroup.50kb.trees -q n82_fossil.wOutgroup.50kb.ASTRAL.spTree -o n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored -t 1 2> n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored.log

# To get the local posterior probability for each branch, score the species tree using option 3:
java -jar $ASTRAL/astral.5.7.8.jar -i n82_fossil.wOutgroup.50kb.trees -q n82_fossil.wOutgroup.50kb.ASTRAL.spTree -o n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored.localPP -t 3 2> n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored.localPP.log

# Rename individuals in newick files before plotting trees, using newich utils:
nw_rename n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored n82_fossil_nw_rename_map.txt > n82_fossil.wOutgroup.50kb.ASTRAL.spTree.scored.rename
## where n82_fossil_nw_rename_map.txt contains old_name new_name (white space separated) on each line.
