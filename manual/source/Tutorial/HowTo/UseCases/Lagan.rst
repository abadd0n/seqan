.. sidebar:: ToC

    .. contents::

.. _how-to-use-cases-lagan:

LAGAN
=====

Learning Objective
 You will learn how to write a simple LAGAN algorithm

Difficulty
  Advanced

Duration
  1h

Prerequisites
  :ref:`tutorial-datastrucures-indices-q-gram-index`, :ref:`tutorial-algorithms-seed-extension`

Background: Alignment of genomic sequences
------------------------------------------
By comparing two sequences of two different species, we can discover new information concerning the conservation of functional units.
These discoveries can assist us us in better understanding and analysing the cause of genetic diseases.
In order the compare two biological sequences we have to calculate a global alignment.
This will show us which operations are necessary to transform one into the other.
Unfortunately it is not advisable to perform a standard global alignment in order to align genomic sequences
for its runtime is determined by the product of the length of both sequences.

Alternatively you can calculate local alignments in order to identify homologous regions.
Doing it this way does not allow us to defer connections between the shared order of the functional units. <======= WUT?!
In order to efficiently calculate a global alignment it is advisable to use heuristics that allow us to identify useful local alignment,
which in turn will be utilized to calculate the global alignment.

One of the most famous examples of that kind of heuristic is the LAGAN-Algorithm.
It consists of three basics steps:
B) Generation of local alignments between the two genomes.
C) Construction of a global map by chaining the identifies segments.
D) Calculation of the optimal alignment within the regions not covered by the local alignments.

The goal of this tutorial is to write a simple version of the LAGAN-Algorithm, which will be extended to work iteratively in the last assignment.
Input will be two FASTA-files containing the genomes and the parameters for the seeding step.
The output will consist of a file containing the alignment.

Building qGram-Index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will be reading two sequences from two different FASTA-files.
At first, our application should create a qGram-Index from the database.
An empty ``qGramIndex`` can simply be created with:

.. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
      :fragment: store

Files can be read from disk with the function :dox:`SeqFileIn#readRecord` that expects a file and two ``StringConcept`` objects.
The contents of different files can be loaded with subsequent calls of ``readRecord``.
As we want the user to specify the files via command line, our application will parse them using the :dox:`ArgumentParser` and store them in an option object.

In your first assignment you need to complete a given code template and implement a function that loads two FASTA-files into two different ``StringConcept`` objects.
We will call them seqH and seqV.
Following this we will create a q-gram-index with variable size based on SeqH.


Assignment 1
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**) and implement a way to read two sequences.
     Use the file paths given in the options object and report an error if the files could not be opened.

     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment1.cpp

   Hint
     * Open STL `std::fstream <http://www.cplusplus.com/reference/iostream/ifstream>`_ objects and use the function :dox:`SeqFileIn#readRecord` .
     * `ifstream::open <http://www.cplusplus.com/reference/iostream/ifstream/open>`_ requires the file path to be given as a C-style string (``const char *``).
     * Use `string::c_str <http://www.cplusplus.com/reference/string/string/c_str>`_ to convert the option strings into C-style strings.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_solution1.cpp
           :fragment: solution



Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**) and implement a way to create a q-gram-index for SeqH while using open addressing.

     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment1.cpp

   Hint
     * use :dox:`OpenAddressingQGramIndex` .
     * use the function :dox:`Shape#resize` .

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_solution1.cpp
           :fragment: solution



		   
		   
		   
		   
		   
		   
		   
		   

Assignments
^^^^^^^^^^

Task 1: Indexaufbau
Im ersten Teil soll eine qGram-Index über die Datenbank aufgebaut werden. Die Größe der q-Gramme wird vom Nutzer übergeben.

Task 2: Seeds finden
Anschließend wird die Query gescanned und die gefundenen Seeds in ein SeedSet eingefügt. Dabei soll zu Beginn die einfache merge Methode verwendet werden.

Task 3: Chaining und Alignment
Aus dem erhaltenen SeedSet soll nun ein globale Chain [5] berechnet werden. Diese wird anschließend in mittels Banded-Chain-Alignment Algorithmus [4] zu einem globalen Alignment zusammengefügt.

Additionally
^^^^^^^^^^^
Option 1: CHAOS-Chaining:
Die Methode zum Zusammenführen von Seeds soll durch das CHAOS-Chaining [6] ersetzt werden, dabei muss auch die Eingabe an das Programm ersetzt werden.

Option 2: Multiple Suchinstanzen
Im originalen LAGAN Algorithmus wird das Alignment wiederholt ausgeführt wobei in jedem Schritt die Parameter für das CHAOS-Chaining verändert werden, sodass die q-Gramme kleiner aber die erlaubte Fehlertoleranz größer wird. Zuerst wird dabei nach langen gut konservierten Regionen gesucht und in jedem weiteren Durchlauf nur noch die Bereiche zwischen den gefundenen Seeds mit relaxierten Parametern gesucht. Die Method soll dahingehen erweitert werden, dass sie ebenfalls mit 3 verschiedene CHAOS-Einstellungen berechnet werden kann. e` to efficiently determine which genes overlap a read alignment.




Extract Gene Intervals
^^^^^^^^^^^^^^^^^^^^^^

Now that the Fragment Store contains the whole annotation tree, we want to traverse the genes and extract the genomic ranges they span.
In the annotation tree, genes are (the only) children of the root node.
To efficiently retrieve the genes that overlap read alignments later, we want to use interval trees, one for each contig.
To construct an interval tree, we first need to collect :dox:`IntervalAndCargo` objects in a string and pass them to :dox:`IntervalTree#createIntervalTree`.
See the interval tree demo in ``demos/interval_tree.cpp`` for more details.
As cargo we use the gene's annotation id to later retrieve all gene specific information.
The strings of ``IntervalAndCargo`` objects should be grouped by ``contigId`` and stored in an (outer) string of strings.
For the sake of simplicity we don't differ between genes on the forward or reverse strand and instead always consider the corresponding intervals on the forward strand.

To define this string of strings of ``IntervalAndCargo`` objects, we first need to determine the types used to represent an annotation.
All annotations are stored in the :dox:`FragmentStore::annotationStore` which is a Fragment Store member and whose type is :dox:`FragmentStore::TAnnotationStore`.
The value type of the annotation store is the class :dox:`AnnotationStoreElement`.
Its member typedefs :dox:`AnnotationStoreElement::TPos` and :dox:`AnnotationStoreElement::TId` define the types it uses to represent a genomic position or the annotation or contig id:

.. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
      :fragment: typedefs

The string of strings of intervals can now be defined as:

.. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
      :fragment: interval

In your second assignment you should use an :dox:`AnnotationTreeIterator AnnotationTree Iterator` to traverse all genes in the annotation tree.
For each gene, determine its genomic range (projected to the forward strand) and add a new ``TInterval`` object to the ``intervals[contigId]`` string, where ``contigId`` is the id of the contig containing that gene.

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more..**).
     Implement the function ``extractGeneIntervals`` that should extract genes from the annotation tree (see :dox:`AnnotationTreeIterator AnnotationTree Iterator`) and create strings of :dox:`IntervalAndCargo` objects - one for each config - that contains the interval on the forward contig strand and the gene's annotation id.

     .. container:: foldable

        Extend the definitions:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment2.cpp
           :fragment: definitions

        Add a function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment2.cpp
           :fragment: yourcode

        Extend the ``main`` function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment2.cpp
           :fragment: main

        and

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment2.cpp
           :fragment: main2

   Hint
     .. container:: foldable

        You can assume that all genes are children of the root node, i.e. create an :dox:`AnnotationTreeIterator AnnotationTree Iterator`, :dox:`AnnotationTreeIterator#goDown go down` to the first gene and :dox:`AnnotationTreeIterator#goRight go right` to visit all other genes.
        Use :dox:`AnnotationTreeIterator#getAnnotation` to access the gene annotation and :dox:`IteratorAssociatedTypesConcept#value` to get the annotation id.

     .. container:: foldable

        Make sure that you append :dox:`IntervalAndCargo` objects, where ``i1`` < ``i2`` holds, as opposed to annotations where ``beginPos`` > ``endPos`` is possible.
        Remember to ensure that ``intervals`` is of appropriate size, e.g. with

        .. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
              :fragment: resize

        Use :dox:`StringConcept#appendValue` to add a new ``TInverval`` object to the inner string, see :dox:`IntervalAndCargo::IntervalAndCargo IntervalAndCargo constructor` for the constructor.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_solution2.cpp
           :fragment: solution

Construct Interval Trees
^^^^^^^^^^^^^^^^^^^^^^^^

With the strings of gene intervals - one for each contig - we now can construct interval trees.
Therefore, we specialize an :dox:`IntervalTree` with the same position and cargo types as used for the :dox:`IntervalAndCargo` objects.
As we need an interval tree for each contig, we instantiate a string of interval trees:

.. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
      :fragment: tree

Your third assignment is to implement a function that constructs the interval trees for all contigs given the string of interval strings.

Assignment 3
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**).
     Implement the function ``constructIntervalTrees`` that uses the interval strings to construct for each contig an interval tree.
     **Optional:** Use OpenMP to parallelize the construction over the contigs, see :dox:`SEQAN_OMP_PRAGMA`.

     .. container:: foldable


        Extend the definitions:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment3.cpp
           :fragment: definitions

        Add a function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment3.cpp
           :fragment: yourcode

        Extend the ``main`` function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment3.cpp
           :fragment: main

        and

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment3.cpp
           :fragment: main2

   Hint
     .. container:: foldable

        First, resize the string of interval trees accordingly:

        .. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
              :fragment: resize_tree

   Hint
     .. container:: foldable

        Use the function :dox:`IntervalTree#createIntervalTree`.

        **Optional:** Construct the trees in parallel over all contigs with an OpenMP parallel for-loop, see `here <http://developers.sun.com/solaris/articles/openmp.html>`_ for more information about OpenMP.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_solution3.cpp
           :fragment: solution

Compute Gene Coverage
^^^^^^^^^^^^^^^^^^^^^

To determine gene expression levels, we first need to compute the read coverage, i.e. the total number of reads overlapping a gene.
Therefore we use a string of counters addressed by the annotation id.

.. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
      :fragment: reads

For each read alignment we want to determine the overlapping genes by conducting a range query via :dox:`IntervalTree#findIntervals` and then increment their counters by 1.
To address the counter of a gene, we use its annotation id stored as cargo in the interval tree.

Read alignments are stored in the :dox:`FragmentStore::alignedReadStore`, a string of :dox:`AlignedReadStoreElement AlignedReadStoreElements` objects.
Their actual type can simply be determined as follows:

.. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
      :fragment: read_alignment_type

Given the :dox:`AlignedReadStoreElement::contigId`, :dox:`AlignedReadStoreElement::beginPos`, and :dox:`AlignedReadStoreElement::endPos` we will retrieve the annotation ids of overlapping genes from the corresponding interval tree.

Your fourth assignment is to implement the count function that performs all the above described steps.
Optionally, use OpenMP to parallelize the counting.

Assignment 4
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**).
     Implement the function ``countReadsPerGene`` that counts for each gene the number of overlapping reads.
     Therefore determine for each :dox:`AlignedReadStoreElement` begin and end positions (on forward strand) of the alignment and increment the ``readsPerGene`` counter for each overlapping gene.

     **Optional:** Use OpenMP to parallelize the function, see :dox:`SEQAN_OMP_PRAGMA`.

     .. container:: foldable

        Extend the definitions:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment4.cpp
           :fragment: definitions

        Add a function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment4.cpp
           :fragment: yourcode

        Extend the ``main`` function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment4.cpp
           :fragment: main

        and

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment4.cpp
           :fragment: main2

   Hint
     .. container:: foldable
        First, resize and zero the string of counters accordingly:

        .. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
              :fragment: resize_reads

        Make sure that you search with :dox:`IntervalTree#findIntervals` where ``query_begin < query_end`` holds, as opposed to read alignments where ``beginPos`` > ``endPos`` is possible.

   Hint
     .. container:: foldable

        The result of a range query is a string of annotation ids given to :dox:`IntervalTree#findIntervals` by-reference:

        .. includefrags:: demos/tutorial/simple_rna_seq/base.cpp
              :fragment: result

        Reuse the result string for multiple queries (of the same thread, use ``private(result)`` for OpenMP).

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_solution4.cpp
           :fragment: solution


Output RPKM Values
^^^^^^^^^^^^^^^^^^

In the final step, we want to output the gene expression levels in a normalized measure.
We therefore use **RPKM** values, i.e. the number of **r**\ eads **p**\ er **k**\ ilobase of exon model per **m**\ illion mapped reads (1).
One advantage of RPKM values is their independence of the sequencing throughput (normalized by total mapped reads), and that they allow to compare the expression of short with long transcripts (normalized by exon length).

The exon length of an mRNA is the sum of lengths of all its exons.
As a gene may have multiple mRNA, we will simply use the maximum of all their exon lengths.

Your final assignment is to output the RPKM value for genes with a read counter ``> 0``.
To compute the exon length of the gene (maximal exon length of all mRNA) use an :dox:`AnnotationTreeIterator AnnotationTree Iterator` and iterate over all mRNA (children of the gene) and all exons (children of mRNA).
For the number of total mapped reads simply use the number of alignments in the :dox:`FragmentStore::alignedReadStore`.
Output the gene names and their RPKM values separated by tabs as follows:

.. includefrags:: demos/tutorial/simple_rna_seq/genequant_solution5.cpp.stdout


.. todo: Move the files to somewhere else.

Download and decompress the attached mouse annotation (`Mus_musculus.NCBIM37.61.gtf.zip <http://ftp.seqan.de/manual_files/seqan-1.4/Mus_musculus.NCBIM37.61.gtf.zip>`_ and the alignment file of RNA-Seq reads aligned to chromosome Y (`sim40mio_onlyY.sam.zip <http://ftp.seqan.de/manual_files/seqan-1.4/sim40mio_onlyY.sam.zip>`_).
Test your program and compare your output with the output above.

Assignment 5
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**).
     Implement the function ``outputGeneCoverage`` that outputs for each expressed gene the gene name and the expression level as RPKM as tab-separated values.

     .. container:: foldable

        Add a function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment5.cpp
           :fragment: yourcode

        Extend the ``main`` function:

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_assignment5.cpp
           :fragment: main

   Hint
     .. container:: foldable

        To compute the maximal exon length use three nested loops: (1) enumerate all genes, (2) enumerate all mRNA of the gene, and (3) enumerate all exons of the mRNA and sum up their lengths.

   Hint
     .. container:: foldable

        Remember that exons are not the only children of mRNA.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/simple_rna_seq/genequant_solution5.cpp
           :fragment: solution

Next Steps
----------

* See :cite:`Mortazavi2008` for further reading.
* Read the :ref:`tutorial-io-sam-bam-io` Tutorial and change your program to stream a SAM file instead of loading it as a whole.
* Change the program such that it attaches the RPKM value as a key-value pair (see :dox:`AnnotationTreeIterator#assignValueByKey`) to the annotation of each gene and output a GFF file.
* Continue with the :ref:`tutorial` rest of the tutorials]].