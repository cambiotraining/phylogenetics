---
title: "Molecular Phylogenetics - an Introduction"
---

::: {.callout-tip}
## Learning Objectives

- Fundamental concepts in phylogenetics.
- What is a sequence alignment and how to estimate one.
- Distance-based phylogenetic methods.
- Maximum likelihood phylogenetics.
- Bayesian phylogenetics.
<!-- - Common pitfalls of phylogenetic inference. -->

:::



## Introduction

The evolutionary history of species is determined in large part by speciations. A phylogenetic tree summarizes such a history, by representing speciations as bifurcations. Each branch of the tree represents evolution of an ancestral species.

Similarly, phylogenetic trees can be used to describe short-term evolutionary histories, such as those of pathogens within epidemics. In this context, bifurcations in the tree are intertwined with transmission events between hosts.

In all these cases, and more, the phylogenetic tree shapes the relatedness of the considered genomes, with genomes being nearer in the tree typically being more similar to each other. Using this fact, we can reconstruct phylogenetic trees starting from observed genomes - with this course, we will discuss how, and how to prevent common pitfalls in molecular phylogenetics.


## What is a phylogenetic tree

A phylogenetic tree is a graph: it is made of nodes and edges, with one edge connecting two nodes.

A node can represent an extant species, and extinct one, or a sampled pathogen: these are all cases of "terminal" nodes, nodes in the tree connected to only one edge, and usually associated with data, such as a genome sequence.

A tree also contains "internal" nodes: these usually represent most recent common ancestors (MRCAs) of groups of terminal nodes, and are typically not associated with observed data, although genome sequences and other features of these ancestors can be statistically inferred. An internal node is most often connected to 3 branches (two descendants and one ancestral), but a multifurcation node can have any number >2 of descendant branches.

::: {.callout-warning}
## Multifurcations
A multifurcation can happen for example when a species radiates, or speciates in >2 descendant species at the same time. Similarly, multifurcations in pathogen trees can represent superspreading events. However, most often multifurcations are inferred from data when a branch of the true evolutionary history tree is not long enough to have caused mutations in the genetic sequences being considered.

When dealing with multifurcations, be aware that many phylogenetic software will not allow them, which might cause multifurcating trees to be represented as bifurcating trees with very short branches, or might cause errors in the execution of some software.
:::

An essential feature of phylogenetic trees is that they have no cycles: for a rooted tree, this means that each non-root node has exactly one parent branch above it connected to it. For an unrooted tree, this means that you cannot find a path starting and ending at the same node and never traversing the same branch twice.

::: {.callout-warning}
## Lateral gene transfers
The assumption of no cycles in a tree means that, in practice, one phylogenetic tree in itself cannot describe hybridizations and other instances of "lateral" evolution such as recombination, gene transfers, etc, and in general all cases in which one genome can have multiple parents.
These events can instead be described with collections of trees (for example one tree for each portion of the genome that shares the same phylogenetic tree, or "gene tree"), or with generalizations of the tree structure, such as phylogenetic networks, ancestral recombination graphs, etc.
:::


### Branch lengths and timed trees

The branches of a tree can be associated with numerical values: their lengths. These lengths however do not always represent the same thing, and in particular one has to be careful with the unit in which branch lengths are expressed.

Most commonly, branch lengths represent divergence. Divergence can be represented in terms of number of mutations per genome position expected on a branch: for example a branch length of 0.01 will mean that, approximately, we expect that 1% of the genome of the considered species/lineage mutated along this branch.
However, some phylogenetic methods, in particular parsimony-based methods, might represent branch lengths with integer numbers, for example a value of 24 would mean 24 substitutions expected on this branch for the considered genetic sequence.

However, branch lengths can also represent time, and trees with such branch lengths are called "timed trees".
While we can estimate divergence-based branches directly from the genome data, to infer timed trees we need additional information.
A species tree can be timed using information about the molecular clock, or about the time of fossils.
A pathogen tree is often timed using the time data of when pathogen samples were collected, since older samples are expected to be closer in time distance to the root of the tree.


### Tree rooting

A phylogenetic tree can be rooted: this means that a root node is specified, and in turn this allows to determine which nodes in the tree are ancestors of which other nodes.
In an unrooted tree this is not possible.
To make a tree rooted ("tree rooting") one often uses:
- An outgroup,
- The assumption that the tree is ultrametric (all terminal nodes have the same distance from the root),
- Time information (such as the time in which the pathogen samples were collected),
- A non-reversible model (genome evolution can leave a directional mark on the genome, for example if the number of T nucleotides increases with time),
- Some or all of the above.


### Tree topology

A clade is the set of all terminal nodes descending from the same ancestor. Each branch and internal node in a tree is associated with a clade. If two trees have the same clades, we say that they have the same topology. If they have the same clades and the same branch lengths, the two tree are equivalent, that is, they represent the same evolutionary history. 


### Newick format

We often need to represent trees in text format, for example to communicate them as input or output of phylogenetic inference software.
The Newick format is the most common text format for phylogenetic trees.

The Newick format encloses each subtree (the part of a tree relating the terminal nodes part of the same clade) with parenthesis, and separates the two child nodes of the same internal node with a ",". At the end of a Newick tree there is always a ";".

For example, the Newick format of a rooted tree relating two samples "S1" and "S2", with distances from the root respectively of 0.1 and 0.2, is

`(S1:0.1,S2:0.2);`

If we add a third sample "S3" as an outgroup, the tree might become

`((S1:0.1,S2:0.2):0.3,S3:0.4);`


### Applications of phylogenetic trees

In may cases, the phylogenetic tree represents the end results of an analysis, for example if we are interested in the evolutionary history of a set of species.

However, in many cases a phylogenetic tree represents an intermediate step, and there are many ways in which phylogenetic trees can be used to help understand evolution and the spread of infectious disease.

In many cases, we may want to know more about genome evolution, for example about mutational pressures, but more frequently about selective pressures. Selection can affect genome evolution in many ways such as slowing down evolution of portion of the genome in which changes are deleterious ("purifying selection"). Instead, "positive selection" can favor changes at certain positions of the genome, effectively accelerating their evolution. Using genome data and phylogenetic trees, molecular evolution methods can infer different types of selection acting in different parts of the genome and different branches of a tree (see e.g. https://doi.org/10.1146/annurev.genet.39.073003.112420 ).

In other scenarios, we might want to combine genetic sequence data with geographical data to reconstruct how species, or pathogens, spread through space and time. This is what we call phylogeography (see e.g. https://doi.org/10.1146/annurev.ecolsys.38.091206.095702 , https://doi.org/10.1016/j.coviro.2011.10.003 ).

It is also possible to use the shape of a phylogenetic tree to reconstruct the past prevalence of an infectious disease - this is called "phylodynamics" (see e.g. https://doi.org/10.1371/journal.pcbi.1002947 ).

This is only a few of many possible applications!


### Caveats

As mentioned before, any non tree-like evolutionary process cannot be captured by phylogenetic trees and has the potential of misleading phylogenetic analyses.

As genes evolve, they can undergo gene loss or duplication. This is one of several reasons why gene trees (where splits might represent gene duplication events and not just speciation events) can differ from species trees.

With hybridation a new species can be formed by merging the genomes of other species. This means that part of the genes will come from one ancestor and the rest from another. Likewise, part of the gene trees will support one species tree, and the rest will support another one.

The genomes of many Eukaryotes undergo recombination each generation. This means that one usually cannot build a single tree from the multiple genomes within the same species; instead, one can infer gene trees for smaller portions of the genome.
When using a single genome from each species, it is common to assume a single tree. However, even in this case within-species recombination of ancestral species can cause diversity in gene trees, a phenomenon we call "incomplete lineage sorting" (see e.g. https://doi.org/10.1093/sysbio/syw082 ).

With lateral gene transfer, genes move from one species to another. These genes will then have different phylogenetic trees.

In bacteria, recombination has usually a more limited effect than in Eukaryotes. As such, while different parts of the genome might have evolved under different trees, usually these trees are similar to the backbone phylogeny of the non-recombining parts of the genome, which we call "clonal frame", see e.g. https://doi.org/10.1016/j.tim.2010.04.002 .




## Alignment

In molecular phylogenetic we use sequence data to infer trees. However, sequence data usually needs to be "aligned". This is because the same homologous sequence can be of different length for different species due to indels, and this makes it hard to reconstruct which substitutions occurred where.
An alignment is a matrix containing one sequence per row. To make the rows/sequences in the matrix of the same length, we pad them with gap "-" characters.

However, we can't put arbitrarily gap characters in the matrix: this would lead to inaccurate trees.
Say we have sequences 
<br />
`ATCG`
`CG`

while this might be a reasonable alignment:
<br />
`ATCG`
<br />
`--CG`

this one is typically much more unlikely:
<br />
`A-TCG`
<br />
`-CG--`

Alignment methods are usually guided by one of two principles:
- Evolutionary alignments want to place homologous characters in the same alignment column. Characters are homologous if they descended from the same ancestral character: while they might have been substituted with different characters, they have not been inserted. This is usually the same assumption made by phylogenetic methods when interpreting an alignment.
- Structural alignments want similar characters to be in the same column, so it makes sense that two characters are in the same column even if they are not homologous, but if they have similar functions within the sequence.

An evolutionary alignment represents an inference regarding the evolutionary sequence of the considered species, in particular regarding insertions and deletions. Phylogenetic inference methods usually don't use insertion and deletion events to infer the tree, but they still assume that characters in the same alignment column are homologous.


### Alignment inference

Pairwise alignment means aligning only two sequences. One popular technique for pairwise alignment is dynamic programming in a form similar to the Needleman-Wunsch algorithm.

Multiple sequence alignment (MSA) means aligning >2 sequences, and is a much more computationally complex process.
A common heuristic strategy is "progressive alignment", where 2 sub-alignments at the time are pairwise aligned, starting from the 2 most closely related sequences, as defined in a given guide tree.

With iterative alignment, one sequence at the time is removed from the alignment, and then re-aligned with the rest, with the aim of improving the alignment.

There are many alignment software available, e.g. Clustal, MAFFT, MUSCLE, PRANK, T-Coffee, ProbAlign. At short divergence these typically give similar results, but with highly divergent sequences alignments can be very different between methods. 
There is no clear consensus on which approach is best, although MAFFT is often considered a great compromise of accuracy/speed; for the inference of selection, there is evidence that evolutionary aligners like PRANK are robust to biases leading to falsely infer positive selection due to overalignment. 
There is no clear consensus regarding methods to "clean" alignments or to assess the confidence of different alignment regions. 


### Alignment practical

Given a set of non-aligned sequences, let's generate one or more alignments (to be completed...).



## Phylogenetic inference - distance-based methods

Phylogenetic relationships impact genomes in that closely related sequences are expected to be more similar to each other.
This is the basic idea behind distance-based methods: the phylogenetic distance between two samples (the sum of the branch lengths separating the two sequences in the tree) is expected to be proportional to the divergence between the considered genomes.

Distance-based methods first calculate the distances between pairs of aligned sequences, and from these they estimate their pairwise divergence. These divergence estimates (which make up the "distance matrix") are then used as constraints to build a phylogenetic tree.

Neighbor-Joining (NJ), one of the most popular distance-based approaches, first creates clades of the most closely related sequences. Then, it updates the distance matrix, replacing the entries for these sequences with distances for the newly created clades. This is repeated iteratively until all clades and branch lengths have been estimated.

Distance-based methods are popular thanks to their simplicity and computational efficiency. In particular, NJ trees are often used as starting trees for other more complex phylogenetic inference methods.


### Neighbor-Joining inference practical

Infer a tree with NJ (to be completed...)


## Phylogenetic inference - Maximum parsimony

Maximum parsimony phylogenetic inference is based on the assumption that a tree that can explain the data (the alignment) requiring the fewest mutations is the preferable tree.

### Tree search

Maximum parsimony methods usually start from an initial tree, which can be estimated for example with NJ, and then try to iteratively improve the tree by attempting to decrease its parsimony score.
The parsimony score is defined as the minimum number of mutations required by a tree to explain the alignment.
One cannot attempt to evaluate all possible trees, because their number is typically too large, and this is why we use heuristic tree searches.

To improve the tree, methods often perform small changes to the tree, and then evaluate the parsimony score of the modified tree.
A type of change is "nearest neighbor interchange", where two nearby subtrees are swapped. Another is "subtree prune and regraft" when a subtree is removed from the tree and re-attached somewhere else.

### Long branch attraction

While maximum parsimony trees are generally regarded as useful, they are not exempt from biases. One of the most studied of these biases is "long branch attraction": when very long branches are present in the true tree, these branch tend to be biasedly clustered together in the inferred maximum parsimony trees. These biases are however not expected when studying short divergence data.

### Maximum parsimony inference practical

Infer a tree with maximum parsimony?




## Phylogenetic inference - maximum likelihood

The phylogenetic likelihood of a tree T is the probability of the alignment A conditional on the tree T: P(A|T).
Maximum likelihood methods use the likelihood function to evaluate how realistic ("likely") different trees are, and aim at returning the tree with the highest value (the maximum likelihood tree).

The likelihood function is typically computed by assuming that all sites in the alignment evolve independently of one another, and by using a dynamic programming approach (the Felsenstein pruning algorithm) to calculate the likelihood of each site.

### Substitution models

A key component need for calculating a tree likelihood is a substitution model. Substitution models describes how states (be it nucleotides, amino acids, codons, or phenotypes) transition between one another. For example, two nucleotides with a high substitution rates between them are expected to substitute each other more frequently.
Substitution models can be broadly split between parametric (the entries of the substitution matrix are estimated from the considered data) or empirical (the rates are pre-estimated, usually using a large database).
Model selection methods can select a substitution model for phylogenetic inference according to different optimality criteria.

### Tree search

Maximum likelihood methods, like maximum parsimony methods, start from an initial tree, and then try to iteratively improve it by attempting to increase its likelihood score typically via "nearest neighbor interchange" and "subtree prune and regraft" moves.

It is not guarantee that these sets of move will lead us to find the maximum tree - we might instead land on a local maximum tree, one that is not optimal, but that we can't improve anymore with our set of moves.


### Maximum likelihood inference practical

Infer a tree with maximum likelihood? Maybe also run model selection software?


## Tree uncertainty - bootstrap

All the methods for phylogenetic inference that we discussed so far aim at estimating a single realistic tree, but they don't automatically tell us how confident we should be in the tree, or in individual branches of the tree.

One common way to address this limitation is using the phylogenetic bootstrap approach (Felsenstein, 1985).
This consist first in sampling a large number (say, 1000) of bootstrap alignments.
Each of these alignments has the same size as the original alignment, and is obtained by sampling with replacement the columns of the original alignment; in each bootstrap alignment some of the columns of the original alignment will usually be absent, and some other columns would be represented multiple times.
We then infer a bootstrap tree from each bootstrap alignment.
Because the bootstrap alignments differ from each other and from the original alignment, the bootstrap trees might different between each other and from the original tree.
The bootstrap support of a branch in the original tree is then defined as the proportion of times in which this branch is present in the bootstrap trees.


## Bayesian phylogenetic inference

Bayesian methods, like maximum likelihood ones, also use the likelihood function to score trees.
However, their aim is not to retrieve the tree with the highest likelihood, but rather to sample from the posterior distribution of trees.

This means that trees are not evaluated solely based on their likelihood, but also based on a "tree prior" function. These tree priors can be very useful for estimating model parameters that affect the tree shape, for example the history of prevalence of an infectious disease through time.

It also means that the output of a Bayesian method includes multiple trees, representing tree inference uncertainty.

Bayesian phylogenetics is particularly useful at short divergence, since parts of the of phylogeny might be very uncertain due to lack of substitutions.

Bayesian methods typically use Monte Carlo Markov Chain (MCMC) for tree space exploration.
This is somewhat similar to the maximum likelihood tree search, in that a current tree is randomly modified, and the modification can be accepted or rejected based on its probability. 
However, MCMC allows the tree probability to decrease as well as to increase. 
At the end of the MCMC run, an initial number of trees is discarded (the "burn-in") and a sample of the trees and parameter values visited is given to the user.


### Bayesian inference practical

Use BEAST 1 or 2 for a simple phylogenetic estimation run?


## Summary

- Molecular phylogenetics typically requires aligned sequences. Several methods are available for sequence alignment.
- Distance-based methods like neighbor-joining offer efficient ad heuristic tree estimates.
- Maximum likelihood methods are typically slower but are statistically consistent and can include complex substitution models.
- Bayesian methods are usually the slowest but most encompassing ones, measuring tree uncertainty and allowing complex evolutionary models such as phylodynamics.
