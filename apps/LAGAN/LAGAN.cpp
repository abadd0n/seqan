//
// Created by abadd0n on 4/30/16.
//

#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>

// fragen:  wie funktioniert der kack
//          how do we find the best seed set?
//
//

using namespace seqan;
template <typename TSeedChain, typename TSeqString>
inline bool createSeedChain(TSeedChain & seedChain, TSeqString const & seqH, TSeqString & seqV, unsigned const & qgramSize)
{
    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;
    typedef Infix<Dna5String>::Type TInfix;

    TIndex index(seqH);
    resize(indexShape(index),qgramSize);
    TInfix kmer;

    //create seedSet
    typedef Seed<Simple>        TSeed;
    typedef SeedSet<TSeed>      TSeedSet;


    TSeedSet seedSet;
    Score<int, Simple> scoringScheme( 2,-1,-2 ); // match, mismatch, gap
    for (unsigned qPos = 0; qPos < length(seqV)-qgramSize+1; ++qPos)
    {
        kmer = infix(seqV, qPos, qPos+qgramSize);
        hash(indexShape(index), begin(kmer));
        for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i) {
            // add seed to seed set using CHAOS chaining method
            if ( !addSeed( seedSet, TSeed(getOccurrences( index, indexShape(index) )[i], qPos, qgramSize), 5 /*max diag dist*/, 10 /*band width*/, scoringScheme, seqH, seqV, Chaos() ) )
                addSeed(seedSet, TSeed(getOccurrences( index, indexShape(index) )[i], qPos, qgramSize), Single()); // just add seed if CHAOS fails

            // Probably we don't need to extend them manually
            //extendSeed( seed, d5s1, d5s2, EXTEND_BOTH, scoringScheme, 3 /*score drop off*/, UngappedXDrop() );
            //appendValue(seedChain, TSeed(getOccurrences(index, indexShape(index))[i], qPos, qgramSize));
        }
    }


    chainSeedsGlobally(seedChain, seedSet, SparseChaining());

    return 1;
}

//int main(int argc, char *argv[]) {
int main()
{
    CharString filePath1 = getAbsolutePath("apps/LAGAN/bovine_adenovirus_6.fa");
    CharString filePath2 = getAbsolutePath("apps/LAGAN/bovine_adenovirus_D.fa");
    SeqFileIn sfi1;
    if ( !open( sfi1, toCString(filePath1) ) ) {
        std::cerr << "ERROR: Could not open the first file" << std::endl;
        return 1;
    }
    SeqFileIn sfi2;
    if ( !open( sfi2, toCString(filePath2) ) ) {
        std::cerr << "ERROR: Could not open the second file" << std::endl;
        return 1;
    }
    CharString meta;
    Dna5String seqH;
    Dna5String seqV;
    readRecord(meta, seqH, sfi1);
    readRecord(meta, seqV, sfi2);

    const unsigned qgramSize = 15;
    String<Seed<Simple>> seedChain;

    createSeedChain(seedChain, seqH, seqV, qgramSize);




    Score<int, Simple> scoringSchemeAnchor(0, -1, -1);
    Score<int, Simple> scoringSchemeGap(2, -1, -1, -2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seqH);
    assignSource(row(alignment, 1), seqV);
    AlignConfig<true, true, true, true> alignConfig;

    int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;


    return 0;
}