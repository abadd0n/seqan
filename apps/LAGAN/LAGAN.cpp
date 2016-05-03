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
//          it doesn't work! WHYY?!
//          IT DOES WORK!! WHYYYY?!?!?!

using namespace seqan;
template <typename TSeedString>
bool printSeedString(TSeedString & seedSet) {
    typedef typename Iterator<TSeedString>::Type            TSeedChainIter;
    std::cout << "SeedSet" << std::endl;
    for (TSeedChainIter it = begin(seedSet); it != end(seedSet); ++it) {
        if ( endPositionH(*it) < beginPositionH(*it) || endPositionV(*it) < beginPositionV(*it) )
            return false;
        std::cout << "SeqH\tStart:\t" << beginPositionH(*it) << "\tEnd:\t" << endPositionH(*it) << "\tSeqV\tStart:\t" << beginPositionV(*it) << "\tEnd:\t" << endPositionV(*it) << std::endl;
    }
    return true;
}

template <typename TSeedChain>
void updateSeedPositions(TSeedChain & seedChain, unsigned globalPosH, unsigned globalPosV)
{
    typedef typename Iterator<TSeedChain>::Type            TSeedChainIter;
    for (TSeedChainIter it = begin(seedChain); it != end(seedChain); ++it) {
        setBeginPositionH(*it, beginPositionH(*it) + globalPosH);
        setEndPositionH(*it, endPositionH(*it) + globalPosH);
        setBeginPositionV(*it, beginPositionV(*it) + globalPosV);
        setEndPositionV(*it, endPositionV(*it) + globalPosV);
    }
}

template <typename TSeedChain, typename TSeqString, typename TQGramSize, typename TSeed>
inline bool createSeedChain(TSeedChain & seedChain, TSeqString const & seqH, TSeqString & seqV, TQGramSize const & qgramSize, TSeed const & /*tag*/)
{
    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<SimpleShape, OpenAddressing > > TIndex;
    typedef Infix<Dna5String>::Type TInfix;

    TIndex index(seqH);
    resize(indexShape(index),qgramSize);
    TInfix kmer;

    //create seedSet
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

    if ( empty (seedChain) )
        return false;

    return true;
}

//int main(int argc, char *argv[]) {
int main() {
    CharString filePath1 = getAbsolutePath("apps/LAGAN/bovine_adenovirus_6.fa");
    CharString filePath2 = getAbsolutePath("apps/LAGAN/bovine_adenovirus_D.fa");
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
    CharString meta;
    Dna5String seqH;
    Dna5String seqV;
    readRecord(meta, seqH, sfi1);
    readRecord(meta, seqV, sfi2);

    const unsigned qgramSizes[] = {10, 3};
    String<Seed<Simple>> seedChain;

    createSeedChain(seedChain, seqH, seqV, qgramSizes[0], Seed<Simple>());

    typedef Iterator<String<Seed<Simple> > >::Type TSeedChainIter;
    String<Seed<Simple>> localSeedChain;
    std::cout << "Length seqH:\t" << length(seqH) << "\tseqV:\t" << length(seqV) << std::endl;

    {
        unsigned i = 0;
        unsigned pos = 0;
        unsigned infixVBegin;
        unsigned infixHBegin;
        unsigned infixVEnd = 0;
        unsigned infixHEnd = 0;
        unsigned nextInfixVBegin = 0;
        unsigned nextInfixHBegin = 0;

        std::cout << "seedChain length:\t" << length(seedChain) << std::endl;

        unsigned const seedChainLen = length(seedChain);

        for (TSeedChainIter it = begin(seedChain); i < seedChainLen+1; ++it)
        {

            // TODO check if global seed alignment is sorted
            // TODO what does the (*it)?
            // TODO can we use index() with iterators instead of seqs?

            if (it == end(seedChain)) {
                break;
            }

            clear(localSeedChain);

            infixHBegin = nextInfixHBegin;
            infixVBegin = nextInfixVBegin;
            infixHEnd = beginPositionH(*it);
            infixVEnd = beginPositionV(*it);
            nextInfixHBegin = endPositionH(*it);
            nextInfixVBegin = endPositionV(*it);

            if (i == seedChainLen) {
                infixHEnd = length(seqH);
                infixVEnd = length(seqV);
            }

            std::cout << "seqH start:\t" << infixHBegin << "\tend:\t" << infixHEnd << std::endl;
            std::cout << "seqV start:\t" << infixVBegin << "\tend:\t" << infixVEnd << std::endl;




            //std::cout << infix( seqV, infix, beginPositionV( *(it+1) ) ) << std::endl;
            Dna5String seqVInfix = infix(seqV, infixVBegin+1, infixVEnd-1);
            Dna5String seqHInfix = infix(seqH, infixHBegin+1, infixHEnd-1);
            if (createSeedChain(localSeedChain,
                                seqHInfix,
                                seqVInfix,
                                qgramSizes[1], Seed<Simple>())) {
                //append(seedChain,localSeedChain);
                updateSeedPositions(localSeedChain, infixHBegin, infixVBegin);
                insert(seedChain, pos, localSeedChain);
                pos += length(localSeedChain);
                it += length(localSeedChain);
            }
            ++i; ++pos;
        }
    }

    if (!printSeedString(seedChain)){
        std::cerr << "end before begin in one of the seeds" << std::endl;
        return 1;
    }


    Score<int, Simple> scoringSchemeAnchor(0, -1, -1);
    Score<int, Simple> scoringSchemeGap(2, -1, -1, -2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seqH);
    assignSource(row(alignment, 1), seqV);
    AlignConfig<false, false, false, false> alignConfig;

    int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;


    return 0;
}