//
// Created by abadd0n on 4/30/16.
//

#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>


using namespace seqan;

int main(int argc, char *argv[]) {
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

    std::cout << d5s1 << std::endl;
    std::cout << d5s2 << std::endl;


    unsigned qgramSize = 3;
    CharString database = "CATGATTACATA";

    //create qgram index with size specified in variable qgramsize
    typedef Index<DnaString, IndexQGram<UngappedShape<qgramSize> > > TIndex;
    TIndex index(d5s1);

    // Default-construct seed.
    Seed<Simple> seed1;
    std::cout << "seed1\n"
    << "beginPositionH == " << beginPositionH(seed1) << "\n"
    << "endPositionH == " << endPositionH(seed1) << "\n"
    << "beginPositionV == " << beginPositionV(seed1) << "\n"
    << "endPositionV == " << endPositionV(seed1) << "\n"
    << "lowerDiagonal == " << lowerDiagonal(seed1) << "\n"
    << "upperDiagonal == " << upperDiagonal(seed1) << "\n\n";

    // Construct seed with begin and end position in both sequences.
    Seed<Simple> seed2(3, 10, 7, 14);
    setUpperDiagonal(seed2, -7);
    setLowerDiagonal(seed2, -9);
    std::cout << "seed2\n"
    << "beginPositionH == " << beginPositionH(seed2) << "\n"
    << "endPositionH == " << endPositionH(seed2) << "\n"
    << "beginPositionV == " << beginPositionV(seed2) << "\n"
    << "endPositionV == " << endPositionV(seed2) << "\n"
    << "lowerDiagonal == " << lowerDiagonal(seed2) << "\n"
    << "upperDiagonal == " << upperDiagonal(seed2) << "\n\n";



    return 0;
}