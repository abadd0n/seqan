//![assignment]
#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <seqan/score.h>
#include <seqan/sequence.h>

using namespace seqan;

//![updateSeedPositions]
template <typename TSeedChain>
void updateSeedPositions(TSeedChain & seedChain, unsigned globalPosH, unsigned globalPosV)
//![updateSeedPositions]
{
    // INSERT YOUR CODE HERE ...
    //
}

//![createSeedChain]
template <typename TSeedChain, typename TSeqString, typename TQGramSize, typename TSeed>
inline bool createSeedChain(TSeedChain & seedChain, TSeqString const & seqH, TSeqString & seqV, TQGramSize const & qgramSize, TSeed const & /*tag*/)
//![createSeedChain]
{
    // INSERT YOUR CODE HERE ...
    //
}

//int main(int argc, char *argv[]) {
int main() {

    // read sequences from fasta files

    // create file objects for the fasta files
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

    // read the first sequences from each file object
    CharString meta;
    Dna5String seqH;
    Dna5String seqV;
    readRecord(meta, seqH, sfi1);
    readRecord(meta, seqV, sfi2);


    // in this string of seeds we will store the seed chains we will find
    String<Seed<Simple>> seedChain;

    // define iterator for string of seeds
    typedef Iterator<String<Seed<Simple> > >::Type TSeedChainIter;

    std::cout << "Length seqH:\t" << length(seqH) << "\tseqV:\t" << length(seqV) << std::endl;

    // we seek for seeds several times. Every next iteration happens in windows where we haven't found any seeds yet.
    // we store those seeds found in repetitive iterations in localSeedChain
    String<Seed<Simple>> localSeedChain;
    for ( unsigned qGramSize : qgramSizes )
    {
        unsigned i = 0;
        unsigned pos = 0;
        unsigned infixVBegin;
        unsigned infixHBegin;
        unsigned infixVEnd = 0;
        unsigned infixHEnd = 0;
        unsigned nextInfixVBegin = 0;
        unsigned nextInfixHBegin = 0;

        unsigned const seedChainLen = length(seedChain);

        for ( TSeedChainIter it = begin( seedChain ); i < seedChainLen + 1; ++it )
        {
            clear(localSeedChain);

            infixHBegin = nextInfixHBegin;
            infixVBegin = nextInfixVBegin;


            if (i == seedChainLen) {
                infixHEnd = length(seqH);
                infixVEnd = length(seqV);
            } else {
                infixHEnd = beginPositionH(seedChain[pos]);
                infixVEnd = beginPositionV(seedChain[pos]);
                nextInfixHBegin = endPositionH(seedChain[pos]);
                nextInfixVBegin = endPositionV(seedChain[pos]);
            }

            std::cout << i << "/" << seedChainLen << "\tInfix SeqV\t start:\t" << infixVBegin << "\tend:\t" << infixVEnd << std::endl << std::endl;

            printSeedString(seedChain);

            Dna5String seqVInfix = infix(seqV, infixVBegin, infixVEnd);
            Dna5String seqHInfix = infix(seqH, infixHBegin, infixHEnd);

            if (createSeedChain(localSeedChain,
                                seqHInfix,
                                seqVInfix,
                                qGramSize,
                                Seed<Simple>() ) )
            {
                updateSeedPositions(localSeedChain, infixHBegin, infixVBegin);
                insert(seedChain, pos, localSeedChain);
                pos += length(localSeedChain);
                it += length(localSeedChain);
            }
            ++i; ++pos;
        }
    }

    Score<int, Simple> scoringSchemeAnchor(0, -1, -1);
    Score<int, Simple> scoringSchemeGap(2, -1, -1, -2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seqH);
    assignSource(row(alignment, 1), seqV);
    AlignConfig<false, false, false, false> alignConfig;  // ordinary global alignment

    int result = bandedChainAlignment(alignment, seedChain, scoringSchemeAnchor, scoringSchemeGap, alignConfig, 2);

    std::cout << "Score: " << result << std::endl;
    std::cout << alignment << std::endl;

    std::cout << toCString(seedChain);

    return 0;
}
//![assignment]