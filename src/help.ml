let index = [
	("build",("Builds Wagner trees�[]. Building multiple trees with 
a randomized addition of terminals allows for the evaluation of
many possible tree topologies and generates a diversity of
trees for subsequent analysis. The arguments of the command
build specify the number of trees to be generated and
the order in which terminals are added during a single tree building
procedure. During tree building, POY5 reports in the Current
Job window of the ncurses interface which of the terminal
addition strategies (e.g. as_is or randomized)
is currently used.","build(20) 
Builds 20 Wagner trees randomizing the order of terminal
addition [Note: because the argument randomized 
is specified by default, it can be omitted].build(trees:20,randomized) 
A more verbose version of the previous example. By default a build
is randomized, but in this case the addition sequence is explicitly
set. For the total number of trees, rather than simply specifying
20, the label trees is used. The verbose version
might be desirable to improve the readability of the script.build(all:30) 
Builds 30 Wagner trees, trying all possible addition positions for
all terminals.build(15,as_is) 
Builds the first Wagner tree using the order of terminals in the first
imported data file and generates the remaining
14 trees using random addition sequences.build(branch_and_bound,trees:5) 
Builds trees using branch and bound method and keeps up to
5 optimal trees in memory.build(constraint:\"cstree.tre\") 
Builds trees using using the set of constraints specified by the
consensus tree file \"cstree.tre\".build(trees:100,optimize(model:max_count:5,branch: 
all_branches)) 
Builds 100 trees and optimizes the likelihood model after every 5
taxa are added to the tree. All branch lengths are optimized after
the addition of each taxon to the tree."));
	("calculate_support",("Calculates the requested support values. POY5 implements support
estimation based on resampling methods (Jackknife�[]
and Bootstrap�[]) and Bremer support�[, ]. The Jackknife and Bootstrap support values are
computed as frequencies of clades recovered in trees specified (see below). 
All the arguments of calculate_support
command are optional and their order is arbitrary. For examples
of scripts implementing support measures see tutorials 5.4,
5.5 and 5.6.","calculate_support(bremer) 
Calculates Bremer support values by performing independent searches
for every node for every tree in memory. This is equivalent to
executing calculate_support(), the default setting.calculate_support(bremer,build(trees:0),swap(trees:2)) 
Calculates Bremer support values by performing swapping on 
each tree in memory for every node and keeping up to two
best trees per search round.calculate_support(bremer,build(of_file:\"new_trees\"), 
swap(tbr,trees:2)) 
Calculates Bremer support values by performing TBR swapping on 
each tree in the file new_trees located in the current
working directory for every node and keeping up to two
best trees per search round. calculate_support(bootstrap) 
Calculates Bootstrap support values under default settings. This
command is equivalent to calculate_support(bootstrap:5,build(trees:1,
randomized),swap(trees:1)).calculate_support(bootstrap:100,build(trees:5),swap(trees:1)) 
Calculates Bootstrap support values performing one random resampling with
replacement, followed by 5 Wagner tree builds (by random addition sequence)
and swapping these trees under the default settings of the command 
swap, and keeping one minimum-cost tree. The procedure
is repeated 100 times.calculate_support(jackknife) 
Calculates Jackknife support values under default settings. This command
is equivalent to calculate_support(jackknife:(resample:5, 
remove:36),build(trees:1,randomized),swap(trees:1)). calculate_support(jackknife:(resample:1000,remove:25), 
build(100),swap(tbr,trees:5)) 
Calculates Jackknife support values randomly removing 25 percent
of the characters, building 100 Wagner trees by random addition
sequence, swapping these trees using tbr, and keeping
up to 5 minimum-cost tree in the final swap per swap (totaling up
to 500 stored trees per replicate). The procedure is repeated 1000
times."));
	("clear_memory",("Frees unused memory. Rarely needed, this is a useful command when the
resources of the computer are limited. The arguments are optional and
their order is arbitrary.","clear_memory(s) 
This command frees memory including all alignment matrices but keeping
unused pool of sequences."));
	("cd",("Changes the working directory of the program. This command is useful
when data files are contained in different directories. It also
eliminates the need to navigate into the working directory before
beginning a POY5 session.","cd(\"/Users/username/docs/poyfiles\") 
Changes the current directory to the directory poyfiles in a Mac OSX 
environment. Filenames with spaces between words need to be escaped, e.g. 
“poy files” should be typed as “poy\  files”. When using a PC, the forward 
slashes should be replaced with backslashes."));
	("echo",("Prints the content of the string argument into a specified type of output.
Several types of output are generated by POY5 which are specified by the
“output class” of arguments (see below). If no output-class arguments are
specified, the command does not generate any output.","echo(\"Building_with_indel_cost_1\",info) 
Prints to the output window in the ncurses interface and to the
standard error in the flat interface the message 
Building_with_indel_cost_1.echo(\"Final_trees\",output:\"trees.txt\") 
Prints the string Final_trees to the file trees.txt.echo(\"Initial_trees\",output) 
Prints the string Initial_trees to the output window in the
ncurses interface, and to the standard output (stdout in the 
flat interface)."));
	("exit",("Exits a POY5 session. This command does not accept any argument.
The exit command is equivalent to quit.","exit() 
Quits the program."));
	("fuse",("Performs Tree Fusing�[] on the trees in memory. 
Tree Fusing can be used to escape local optima
by exchanging clades with identical composition of terminals, differing in 
arrangement between pairs of trees. Only one exchange between 
pairs of trees is evaluated during a single iteration.","fuse(iterations:10,replace:best,keep:100,swap()) 
This command executes the following sequence of operations. In the
first iteration, clades of the same composition of terminals are exchanged
between two trees from the pool of the trees in memory. The cost of the
resulting trees is compared to that of the trees in memory and a subset of
the trees containing up to 100 trees of best cost is retained in memory.
During each iteration of fuse, the trees are subjected to swapping under 
the default settings of swap. The entire procedure is 
repeated nine more times.fuse(optimize:(model:never,branch:join_region)) 
This command performs tree fusing, specifying that the likelihood model is 
never optimized after each round of fusing , but that a maximum of five 
branches are optimized each round.fuse(swap(constraint)) 
This command performs tree fusing 
with modified settings for swapping that follows each iteration. Once
a given iteration is completed, a consensus tree of the files in memory
is computed and used as constraint file for subsequent rounds of swapping (see
the argument constraint of the command
swap)."));
	("help",("Reports the requested contents of the help file on screen.","help(swap) 
Prints the description of the command
swap in the POY Output window of the ncurses
interface or to the standard error in the flat interface.help(\"log\") 
Finds every command with text containing the substring log and
prints them in the POY Output window of the ncurses
interface or to the standard error in the flat interface."));
	("inspect",("Retrieves the description of a POY5 file produced by the command save . 
If the description were not specified by the user, inspect reports 
that the description is not available. If the file is not a proper POY5 file format, a 
message is printed in the POY Output window of the ncurses 
interface or to the standard error of the flat interface.","inspect(\"initial_search.poy\") 
Prints the description of the POY5 file initial_search.poy
located in the current working directory in the POY Output
window of the ncurses interface or to the standard error in the flat
interface. If, for example, the file was saved using
the command save (\"initial_search.poy\",\"Results_of_Total_Analysis\"), 
then the output message is: Results_of_Total_Analysis."));
	("load",("Imports and inputs POY5 files created by the command
save. The name of the file to be loaded
is included in the string argument. All the information of
the current POY5 session will be replaced with the contents
of the POY5 file. If the file is not in proper POY5 file
format, an error message is printed in the POY Output
window of the ncurses interface, or the standard error 
in the flat interface. See the description of the command 
save on the POY5 file and its usage.","load(\"initial_search.poy\") 
Reads and imports the contents of the POY5 file
initial_search.poy, located in the current working
directory."));
	("perturb",("Performs branch swapping on the trees currently in memory using 
temporarily modified (“perturbed”) characters. Once a local optimum is
found for the perturbed characters, a new round of swapping using the
original (non-modified) characters is performed. Subsequently, the costs
of the initial and final trees are compared and the best trees are
selected. If there are n trees in memory prior to searching using
perturb, then the n best trees are selected at the end.
For example, if there are 20 trees currently in memory, 20 individual
perturb procedures will be performed (each procedure
starting with one of the 20 initial trees), and 20 final trees are
produced.","perturb(resample:50,iterations:10) 
Performs 10 successive repetitions of random resampling of 50
characters with replacement. Branch swapping is performed
using alternating SPR and TBR, and and keeping one minimum-cost tree
(the default of swap).perturb(iterations:20,ratchet:(0.18,3)) 
Performs 20 successive repetitions of a variant of the ratchet (see
above) by randomly selecting 18 percent of the characters (sequence
fragments) and upweighting them by a factor of 3. Branch swapping is
performed using alternating SPR and TBR, and keeping one
optimal tree (the default of swap).perturb(iterations:1,transform(tcm:(4,3))) 
Transforms the cost regime of all applicable characters to the new 
cost regime specified by transform (cost of substitution 
4 and cost of indel 3). Subsequently a single round of branch swapping is
performed using alternating SPR and TBR, and and keeping one
optimal tree (the default of swap).perturb(ratchet:(0.2,5),iterations:25,swap(tbr,trees:5)) 
Performs 25 successive repetitions of a variant of the ratchet (see
above) by randomly selecting 20 percent of the characters (sequence
fragments) and upweighting them by a factor of 5. Branch swapping is
performed using TBR and keeping up to 5 optimal trees in each iteration.perturb(transform(static_approx),ratchet:(0.2,5),iterations:
25,swap(tbr,trees:5)) 
Transforms all applicable (i.e. dynamic homology sequence
characters) using transform into static characters.
Therefore, the subsequent ratchet is performed at the level of
individual nucleotides (as in the original implementation),
not sequence fragments. Thus, ratchet is performed by
selecting 20 percent of the characters (individual nucleotides) and
upweighting them by a factor of 5. Branch swapping is performed
using TBR and keeping up to 5 optimal trees in each iteration as in
the example above."));
	("pwd",("Prints the current working directory in the POY Output window of
the ncurses interface and the standard error (stderr) of the flat interface.
The command pwd does not have arguments. The default
working directory is the shell’s directory when POY5 started.","pwd() 
This command generates the following message: “The current
working directory is /Users/username/datafiles/”. The actual reported
directory will vary depending on the directory of the shell when
POY5 started, or if it has been changed using the command
cd()."));
	("quit",("Exits POY5 session. This command does not have any arguments
quit is equivalent to the command exit.","quit() 
Quits the program."));
	("read",("Imports data files and tree files. Supported formats include ASN1, Clustal, FASTA,
GBSeq, Genbank, Hennig86, Newick, NewSeq, Nexus, PHYLIP, POY3,
TinySeq, and XML. Filenames must be enclosed in quotes and, if multiple
filenames are specified, they must be separated by commas. All filenames 
read into POY5 must include the appropriate suffix (e.g. .aln, .fas, 
.fasta, .ss, .tre). The exclusion of these suffixes will result in an error such as 
\"Sys_error (\"No such file or directory\")\". The filename must match exactly.","read(\"/Users/andres/data/test.txt\") 
Reads the file test.txt located in the path
/Users/andres/data/.read(\"28s.fas\",\"initial_trees.txt\") 
Reads the file 28s.fas and loads the trees in parenthetical notation
of the file initial_trees.txt.read(\"SSU*\",\"*.txt\") 
Reads all the files with names starting with SSU, and all the
files with the extension .txt. The types of the data files are determined
automatically.read(nucleotides:(\"chel.FASTA\",\"chel2.FASTA\")) 
Reads the files chel.FASTA and chel2.FASTA, containing nucleotide
sequences.read(aminoacids:(\"a.FASTA\",\"b.FASTA\",\"c.FASTA\")) 
Reads the amino acid sequence files a.FASTA, b.FASTA, and
c.FASTA.read(\"hennig1.ss\",\"chel2.FASTA\",aminoacids:(\"a.FASTA\")) 
Reads the Hennig86 file hennig1.ss, the FASTA file chel2.FASTA
containing nucleotide sequences, and the amino acid sequence file 
a.FASTA.
script
read(annotated:(\"ch1.txt\",\"ch2.txt\"),chromosome:(\"ch3.txt\")) 
Reads three files containing chromosome-type sequence data. The sequences in 
two files, ch1.txt and ch2.txt, contain pipes (“  ”) 
separating individual loci, whereas the sequences in the third (ch3.txt), 
are without predefined boundaries. [Note: see tutorial 5.10, which illustrates
the transformation of these two data types in the same analysis.]read(genome:(\"mt_genomes\",\"nu_genomes\")) 
Reads two files containing genomic (multi-chromosomal) sequence data.read(breakinv:(\"BI_run1\",tcm:\"alphabet\",level:1,tiebreaker: 
first) 
Reads the file, BI_run1, which contains data in custom alphabet format. 
This data is defined in the tcm file, alphabet. The heuristic level of 
median sequence calculation is set to 1. If ties among median states are encountered, 
the first will be chosen. read(custom_alphabet:(\"CA_run1\",tcm:\"alphabet\",level:2,
tiebreaker:last)) 
Reads the file, CA_run1, which contains data in custom alphabet format. 
This data is defined in the tcm file, alphabet. The heuristic level of 
median sequence calculation is set to 2. If ties among median states are encountered, 
the last will be chosen. read(prealigned:(\"18s.aln\",tcm:(1,2))) 
Reads the prealigned data file 18s.aln which was generated from the nucleotide file 18s.FASTA
using the the transformation costs 1 for substitutions and 2 for indels.read(prealigned:(nucleotides:(\"*.nex\"),tcm:\"matrix1\")) 
Reads character data from all the Nexus files as prealigned data using the the transformation cost
matrix from the file matrix1."));
	("recover",("Recovers the best trees found during swapping, even if the swap were
cancelled. This command functions only if the argument recover 
were included in a previously executed 
(in the current POY5 session) command swap. Otherwise, it has no effect.","recover() 
If the command swap (executed
earlier in the current POY5 session) contained the argument recover,
for example, swap(tbr,recover), this command will restore the best
trees recovered during swapping."));
	("rediagnose",("Performs a re-optimization of the trees currently in memory. This
function is useful for sanity checks of the consistency of the data.
Its main usage is for the POY5 developers. This command does not have
arguments.","rediagnose() 
See the description of the command."));
	("redraw",("Redraws the screen of the terminal. This command is only used in the ncurses
interface, other interfaces ignore it. redraw clears the
contents of the Interactive Console window but retains the contents
of the other windows. It does not affect the state of the search and the data
currently in memory.","redraw() 
See the description of the command."));
	("rename",("Replaces the name(s) of specified item(s) (characters or terminals). This command allows 
for substituting taxon names and helps merging multiple datasets without modifying the original
data files. More specifically, it can be used, for example, (1) for housekeeping purposes,
when it is desirable to maintain long verbose taxon names (such as catalog or GenBank
accession numbers) associated with the original data files but avoid reporting these 
names on the trees (although see the note on the usage of a “$” in the taxon name below); 
(2) to provide a single name for a terminal in cases where the corresponding
data are stored in different files under different terminal names; and (3) to change an
outdated or invalid terminal name.","rename(terminals,\"synfile\") 
This command renames terminal names
contained in the synonymy file synfile in all subsequently imported data files.rename(terminals,(\"Mytilus_sp\",\"Mytilus_edulis\")) 
This command
renames the terminal taxon Mytilus_sp as Mytilus_edulis in all subsequently
imported data files.rename(terminals,(\"Chamaeleo_weidersheimi\",\"Trioceros_wieder
sheimi\"),(\"Chamaeleo_wiedersheimi\",\"Trioceros_wiedersheimi\")) 
This command renames the terminal taxa Chamaeleo_wiedersheimi 
and Chamaeleo_weidersheimi (a misspelling of the previous name) as 
Trioceros_wiedersheimi in all subsequently imported data files."));
	("report",("Outputs the results of current analysis or loaded data in the POY Output
window of the ncurses interface, the standard output of the flat
interface, or to a file. To redirect the output to a file, the file name in 
quotes and followed by a comma must be included in the argument list
of report. All arguments for report are
optional. This command allows the user to output information concerning the 
characters and terminals, diagnosis, export static homology data, implied 
alignments, trees, as well as other miscellaneous arguments.","SEC END"));
	("run",("Runs POY5 script file(s). The filenames must be included in
quotes and, if multiple files are included, they must be separated by commas.
The script-containing files are executed in the order in which they are listed
in the string argument.
Executing scripts using run is useful in cases when
operations take a long time or many scripts need to be executed automatically,
e.g. when conducting a sensitivity analysis [].
An in depth description of creating and running scripts is provided in 
the Quickstart. There are no default settings of run.","run(\"script1\",\"script2\") 
This command executes POY5 command scripts contained in the files script1
and script2 in the same order as they are listed in the list of arguments of
run. Recall: If the last line of script1 ends in quit or 
exit, POY5 will finish before script2 can be run."));
	("save",("Saves the current POY5 state of the program to a file (POY5 file). The
first, obligatory string argument specified the name of the POY5 file.
The second, optional string argument specifies a string included in the
POY5 file, that can be retrieved using the command inspect .","save(\"alldata.poy\") 
This command stores all the memory contents of the program in the
file alldata.poy located in the current working directory,
as printed by pwd().save(\"alldata.poy\",\"Total_evidence_data\") 
This command performs the same operation as described in the example above,
but, in addition, it includes the string Total_evidence_data 
with the file alldata.poy,
which can later be retrieved using the command inspect .save(\"/Users/andres/alldata.poy\",\"Total_evidence_data\") 
This command performs the same operation as the command described above
with the important difference that the file alldata.poy generated in the
directory /Users/andres/ instead of the current working directory."));
	("search",("search implements a default search strategy that
includes tree building, swapping using TBR, perturbation using
ratchet, and tree fusing. The strategy involves specifying targets for 
a driven search, such as maximum and minimum execution times, 
maximum allowed memory consumption for tree storage, minimum number of times the
shortest tree is found, and an expected cost for the shortest tree. When executing
search using parallel processing, trees are exchanged upon the
completion of the command (after fusing). Because the lowest cost unique trees 
generated are selected and stored at the end of a search 
(defined by the user with max_time), aggressive use of this 
command in a parallel environment consists of including few sequential
search commands that will allow the processes to
exchange trees and add the pool of selected best trees to subsequent 
iterations of the command (see the example for parallel processing).","search(hits:100,target_cost:385,max_time:1:12:13) 
This command will attempt as many builds, swaps, ratchets, and tree
fusings as possible within the specified time of 1 day, 12 hours, and 13
minutes, finding at least 100 hits (whichever occurs first, the time
limit or the number of hits), knowing that the expected cost of the
best hits is at most 385 steps.For Parallel Implementation of search 
search(max_time:0:6:0)
select()
search(max_time:0:6:0)
select()
search(max_time:0:4:0)
select()
This series of commands will attempt as many builds, swaps, ratchets, and tree
fusings as possible within the specified total time of 16 hours. Trees are exchanged among
processors at the end of each search and the best unique trees
are then selected and included in the following
search command.search(max_time:00:48:00,constraint:\"best_tree.tre\") 
This command will attempt as many builds, swaps, ratchets, and tree
fusings within the specified time period of 48 hours. In this example, however, these 
operations are constrained by the tree specified in the file best_tree.tre. search(max_time:00:48:00,visited:\"visited.txt\") 
This command will attempt as many builds, swaps, ratchets, and tree
fusings as possible within the specified time of 2 days. During this 
time, every visited tree and its cost during the local search will be stored in the file
visited.txt."));
	("select",("Specifies a subset of terminals, characters, or trees from those
currently loaded in memory, to use in subsequent analysis.","select(terminals,names:(\"t1\",\"t2\",\"t3\",\"t4\",\"t5\"), 
characters,names:(\"chel.aln:0\")) 
This command selects only terminals t1, t2,
t3, t4, and t5 and use data only from the
fragment 0 contained in the file chel.aln.select(terminals,files:(\"STL_terminals.txt\")) 
This command selects only the terminals specified in the 
file STL_terminals.txt\". In the data files that are 
subsequently imported, taxa that do not appear in this 
terminals file will be excluded from the analysis.select(terminals,missing:30) 
This command excludes from subsequent analyses all the terminals that
have fewer than 30 percent of characters missing. The list of included and excluded
terminals is automatically reported on screen.
scr
select(optimal) 
Selects all optimal (best cost) trees and discards suboptimal trees from
memory. The pool of optimal trees might contain duplicate trees (that can
be removed using unique).select(unique,within:2.0) 
This command selects all topologically unique optimal and suboptimal trees
the cost of which does not exceed that of the best current cost by more than
2. For example, if the best current cost is 49, all unique trees that fall within
the cost range 49–51 are selected."));
	("set",("Changes the settings of POY5. This command performs diverse auxiliary 
functions, from setting the seed of the random number generator, to
selecting a terminal for rooting output trees, to defining character sets
for different partitions.","set(history:10000,seed:45,log:\"mylog.txt\") 
This command increases the size of the history in the ncurses
interface to 10,000 lines, sets the seed of the random number generator to 45,
and initiates a log file mylog.txt, located in the current
working directory.set(root:\"Mytilus_edulis\") 
This commands selects the terminal Mytilus_edulis as the root
for the output trees.set(iterative:exact) 
Turns on the iterative exact algorithm in all the nucleotide
sequence characters. The program will iterate on each vertex of the
tree until no further tree cost improvements can be made.set(iterative:approximate:2) 
Turns on the iterative approximate algorithm in all the nucleotide
sequence characters. The program will iterate either two times, or
until no further tree cost improvements can be made, whichever
happens first.set(iterative:exact:2) 
Same as the previous, but using the exact algorithm.set(codon_partition:(\"coleop\",names:(\"coleoptera_nd2.fasta\"))) 
Sets codon partitioning of the data file coleoptera_nd2.fasta. Three
sets are partitions named coleop1, coleop2, and
coleop3 will be created.set(opt:exhaustive:3) 
Set floating point optimization to a tolerance of 1e-6, and specify that
a maximum of three optimization iterations occur.set(opt:coarse:10) 
Set floating point optimization to a tolerance of 1e-3, and specify that
a maximum of 10 optimization iterations occur."));
	("store",("Stores the current state of POY5 session in memory. The stored information
includes character data, trees, selections, everything. Specifying the
name of the stored state of the search (using the string argument) does
not, however, generate a file under this name that can be examined;
the name is used only to recover the stored state using the command 
use.","store(\"initial_tcm\") 
transform(tcm:(1,1)) 
use(\"initial_tcm\") 
The first command, store, stores the current
characters and trees under the
name initial_tcm. The second command,
transform, changes the cost regime of molecular characters,
effectively changing the data being analyzed. However, the third
command, use, recovers the initial state stored under the
name initial_tcm."));
	("swap",("swap is the basic local search function in POY5. This
command implements a family of algorithms collectively known in systematics as
branch swapping and in combinatorial optimization as hill climbing. 
They proceed by clipping parts of a given tree and
attaching them in different positions. It can be
used to perform a local search from a set of trees loaded in memory.","swap() 
This command performs swapping under default settings.swap(trees:5) 
Submits current trees to a round of SPR followed by TBR. It keeps
up to 5 minimum cost trees for each starting tree.swap(transform(all,(static_approx))) 
Submits current trees to a round of SPR followed by TBR, using
static approximations for all sequence characters.swap(trees:4,transform(all,(static_approx))) 
Submits current trees to a round of SPR followed by TBR, using
static approximations for all characters, keeping up to 4 minimum
cost trees for each starting tree.swap(constraint:(depth:4)) 
Calculates a consensus tree of the files in memory and uses it as
constraint file, then joins at a distance of at most 4 from the breaking
branch. This is equivalent to swap(constraint:(4)).swap(constraint:(file:\"bleh\")) 
Reads the tree in file bleh and use it as constraint for the
search. This is equivalent to swap(constraint:(\"bleh\")). This 
presumes that the file bleh is located in the current working directory.	swap(drifting:(0.5,2.0)) 
Defines the direction of search via drifting, such that there is a 50%% probability of 
replacing the current tree with a new tree of equal cost. For suboptimal trees with 
a cost d greater than the current best tree, the probability of accepting this tree 
is 1 / (2.0 + d). For example the
probability of keeping a suboptimal tree of 3 steps longer = 1 / (2.0 + 3) = 0.2.swap(sectorial:4) 
Submits current trees to a round of SPR followed by TBR. Join will take place
at a distance equal or less than the value of the argument from the broken edge, 
where the distance is the number of edges in the path connecting them. If no 
argument is given, then no distance limit is set.swap(spr,all,optimize:(model:never,branch:join_region)) 
Submits the current trees to a round of SPR swapping. Following each round of 
spr, the model is never optimized, but a maximum of five branches 
(the new edge, and the two edges on either side of the join site) are optimized.swap(recover,timedprint:(5,\"timedprint.txt\")) 
Saves the current best tree to the file timedprint.txt every 5 seconds."));
	("transform",("Transforms the properties of the imported characters from one type into
another. This includes changing costs for indels and substitution,
modifying character weights, converting dynamic into static homology characters,
transforming nucleotide into chromosomal characters, and specifying characters as either 
static or dynamic likelihood characters, among other operations.","SEC END"));
	("use",("Restores from memory the state of a POY5 session (that includes
character data, selections, trees, all other data and specifications)
that had previously been saved during the session using the
command�store . The recalled session replaces the current
session. The string argument specifies the name of the stored state.","store(\"initial_tcm\") 
transform(tcm:(1,1)) 
use(\"initial_tcm\") 
The first command, store, stores the current characters
and trees under the name initial_tcm. The second command,
transform, changes the cost regime of molecular
characters, effectively changing the data being analyzed. However,
the third command, use, recovers the initial state
stored under the name initial_tcm."));
	("version",("Reports the POY5 version number in the output window of the
ncurses interface, or to the standard error in the flat
interface.","version ()"));
	("wipe",("Deletes the data stored in memory (all character data, trees, etc.).","wipe ()"));
]