//
// Created by abadd0n on 4/30/16.
//

#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>


using namespace seqan;

int main() {
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

    return 0;
}