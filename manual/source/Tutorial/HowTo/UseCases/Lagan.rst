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
These discoveries can assist us in better understanding and analysing the cause of genetic diseases.
In order the compare two biological sequences we have to calculate a global alignment.
This will show us which operations are necessary to transform one into the other.
Unfortunately it is not advisable to perform a standard global alignment in order to align genomic sequences,
for its runtime is determined by the product of the length of both sequences.

Alternatively you can calculate local alignments in order to identify homologous regions.
Doing it this way does not allow us to defer connections between the shared order of the functional units. <======= WUT?!
In order to efficiently calculate a global alignment it is advisable to use heuristics that allow us to identify useful local alignment,
which in turn will be utilized to calculate the global alignment.

One of the most famous examples of that kind of heuristic is the LAGAN-Algorithm.

.. figure:: lagan.jpg
   :width: 60%

   LAGAN example

It consists of three basics steps:
B) Generation of local alignments between the two genomes.
C) Construction of a global map by chaining the identified segments.
D) Calculation of the optimal alignment within the regions not covered by the local alignments.

The goal of this tutorial is to write a simple version of the LAGAN-Algorithm, which will be extended to work iteratively in the last assignment.
Input will be two FASTA-files containing the genomes and the parameters for the seeding step.
The output will consist of a standard output containing the alignment.

Building q-gram-index
^^^^^^^^^^^^^^^^^^^^^

We will be reading two sequences from two different FASTA-files.
At first, our application should create a q-gram-index from the database.

Files can be read from disk with the function :dox:`SeqFileIn#readRecord` that expects a file and two ``StringConcept`` objects.
The contents of different files can be loaded with subsequent calls of :dox:`SeqFileIn#readRecord`.
As we want the user to specify the files via command line, our application will parse them using the :dox:`ArgumentParser` and store them in an option object.

.. includefrags:: demos/tutorial/lagan/assignment1.cpp
    :fragment: include

.. includefrags:: demos/tutorial/lagan/assignment1.cpp
    :fragment: sequences

We will call the reference seqH and the query seqV.

In your first assignment we will begin writing a function that will allow us to to create a seed chain.
The function will look something like this:

.. includefrags:: demos/tutorial/lagan/assignment1.cpp
    :fragment: createSeedChainHead

For now you will only need to complete a given code template and implement a way to create a q-gram-index
with variable size based on the reference.
We will use Open Addressing in order to be able to have the user enter a maximum q-gram-size of 31.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**) and implement a way to create a q-gram-index for SeqH while using open addressing.


     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/assignment1.cpp
            :fragment: createSeedChain

   Hint
     .. container:: foldable

       * use :dox:`OpenAddressingQGramIndex`.
       * use the function :dox:`Shape#resize`.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/solution1.cpp
           :fragment: createSeedChain



Creating a Seed-set and a Global-Chain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that we have a q-gram-index we can begin to find our seeds based on the k-meres from our query.
For this we will use an infix with the specified q-gram-size based on SeqV.
We can find the position of an infix in the reference by hashing it.
This will allow us to add the found seeds to a seed-set using the chaos-chaining method.
Adding a seed via chaos-chaining requires a :dox:`Score`.
Examples on how to implement chaos-chaining can be found here: :dox:`SeedSet#addSeed`.

If a seed cannot be added using the chaos-chaining method,
you should add the seed using the simple-merge method in order to create a new anchor.

An empty ``seedSet`` can simply be created with:

.. includefrags:: demos/tutorial/lagan/solution2.cpp
    :fragment: seedSet

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objectives
     * Use the code template below (click **more...**) and complement the previous function so that it creates a seed-set for SeqV.
     * Additionally combine the seed set into a global chain.

     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/assignment2.cpp
           :fragment: createSeedChain

   Hint 1
     .. container:: foldable

       * use the function :dox:`SegmentableConcept#infix`.
       * use the function :dox:`Shape#hash`.
       * use the function :dox:`SeedSet#addSeed`.
       * use :dox:`Score`.

   Hint 2
     .. container:: foldable

       * use the function :dox:`chainSeedsGlobally`.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/solution2.cpp
           :fragment: createSeedChain



Chaining and Alignment
^^^^^^^^^^^^^^^^^^^^^^

Now that we have successfully created a global chain, it is time to defer the global alignment.

This seedChain will the be extended to a global alignment by using the banded-chain-alignment algorithm.
For this we need to :dox:`Align` both sequences by creating an alignment object and specifying the scoringSchemes.

Assignment 3
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**) and implement a way to create an alignment based on the seedSet.


     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/assignment3.cpp
            :fragment: main

   Hint
     .. container:: foldable

       * use the function :dox:`bandedChainAlignment`.
       * use :dox:`Score`.
       * use :dox:`Align`

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/solution3.cpp
            :fragment: solution


Iterative LAGAN
^^^^^^^^^^^^^^^

.. warning::
    This is not trivial.

The last assignment is more complex, because it will cover how to include iterative steps in the LAGAN-algorithm.
This will allow us to increase the accuracy by trying to find new seed-chains within areas previously not covered by the global chain.

The first step is to write a create seedChain function which will take the following parameters:

.. includefrags:: demos/tutorial/lagan/base.cpp
    :fragment: createSeedChain

Additionally a function is needed that allows the correction of the position of the seedChain found in the gap of the global seedChain,
in order to allow us to add the local seedChain to the global one.
Because we will need the end Positions of the last global seed prior to the gap, the function will look something like this.

.. includefrags:: demos/tutorial/lagan/base.cpp
    :fragment: updateSeedPositions


Assignment 4
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code template below (click **more...**) and implement a way to execute LAGAN iteratively.


     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/base_assignment.cpp
            :fragment: assignment

   Hint
     .. container:: foldable

       * it is important to remember where the local seedChain has to be inserted into the global one.
       * .

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/lagan/base.cpp
            :fragment: solution

Next Steps
----------

* Continue with the rest of the :ref:`tutorial`.