//
// Created by abadd0n on 4/30/16.
//

#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>


using namespace seqan;

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
    Dna5String d5s1;
    Dna5String d5s2;
    readRecord(meta, d5s1, sfi1);
    readRecord(meta, d5s2, sfi2);

    //std::cout << d5s1 << std::endl;
    //std::cout << d5s2 << std::endl;


    const unsigned qgramSize = 3;
    Dna5String database = "CATGATTACATA";

    //create qgram index with size specified in variable qgramsize
    typedef Index<Dna5String, IndexQGram<UngappedShape<qgramSize> > > TIndex;
    typedef Infix<Dna5String>::Type TInfix;

    TIndex index(d5s1);
    TInfix kmere;


    //create seedSet
    typedef Seed<Simple> TSeed;

    Score<int, Simple> scoringSchemeAnchor(0, -1, -1);
    Score<int, Simple> scoringSchemeGap(2, -1, -1, -2);

    String<TSeed> seedChain;

    for (unsigned qPos = 0; qPos < length(d5s2)-qgramSize+1; ++qPos)
    {
        kmere = infix(d5s2, qPos, qPos+qgramSize);
        //std::cout << kmere << std::endl;
        hash(indexShape(index), begin(kmere));
        for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
            appendValue(seedChain, TSeed(getOccurrences(index, indexShape(index))[i], qPos, qgramSize));
            //std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;
    }

    //appendValue(seedChain, TSeed(6, 9, 9, 12));
    //appendValue(seedChain, TSeed(11, 14, 17, 16));

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), d5s1);
    assignSource(row(alignment, 1), d5s2);
    AlignConfig<true, false, false, true> alignConfig;

    int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;

    return 0;
}