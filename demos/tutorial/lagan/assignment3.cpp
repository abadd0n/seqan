#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>


//![main]
int main() {
    // read sequences from fasta files

    // create file objects for the fasta files
    CharString filePath1 = getAbsolutePath("example_reference.fa");
    CharString filePath2 = getAbsolutePath("example_query.fa");
    SeqFileIn sfi1;
    if (!open(sfi1, toCString(filePath1))) {
        std::cerr << "ERROR: Could not open the first file" << std::endl;
        return 1;
    }
    SeqFileIn sfi2;
    if (!open(sfi2, toCString(filePath2))) {
        std::cerr << "ERROR: Could not open the second file" << std::endl;
        return 1;
    }

    // read the first sequences from each file object
    CharString meta;
    Dna5String seqH;
    Dna5String seqV;
    readRecord(meta, seqH, sfi1);
    readRecord(meta, seqV, sfi2);

    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;

    TIndex index(seqH);
    resize(indexShape(index),qgramSize);

    //![seedSet]
    //create seedSet
    typedef Seed<Simple>        TSeed;
    typedef SeedSet<TSeed>      TSeedSet;
    TSeedSet seedSet;
    //![seedSet]

    //define infix
    typedef Infix<Dna5String>::Type TInfix;
    TInfix kmer;

    //scoring scheme defined for chaos-chaining
    Score<int, Simple> scoringScheme( 2,-1,-2 ); // match, mismatch, gap
    for (unsigned qPos = 0; qPos < length(seqV)-qgramSize+1; ++qPos)
    {
        kmer = infix(seqV, qPos, qPos+qgramSize);
        hash(indexShape(index), begin(kmer));
        for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i) {
            // add seed to seed set using CHAOS chaining method
            if ( !addSeed( seedSet,
                           TSeed(getOccurrences( index, indexShape(index) )[i],qPos, qgramSize),
                           5 /*max diag dist*/,
                           10 /*band width*/,
                           scoringScheme,
                           seqH,
                           seqV,
                           Chaos() ) )
                addSeed(seedSet,
                        TSeed(getOccurrences( index, indexShape(index) )[i], qPos, qgramSize),
                        Single()); // just add seed if CHAOS fails
        }
    }

    // INSERT YOUR CODE HERE ...
    //

    return 0;
}
//![main]