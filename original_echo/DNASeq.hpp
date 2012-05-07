#ifndef __DNASEQ_H__
#define __DNASEQ_H__

#include <string>
#include <stdexcept>
#include <tr1/unordered_map>
#include <tr1/memory>
#include <algorithm>

class DNASeq {
    friend class Kmer;
    friend class Kmer_cmp;
    friend class Kmer_hash; 

    static std::tr1::unordered_map<char, char> ComplementBase;

    const char* seq;
    int seq_len;

    static std::tr1::unordered_map<char, char> ComplementBaseInitializer() {
        std::tr1::unordered_map<char, char> ComplementBase;
        ComplementBase['A'] = 'T';
        ComplementBase['C'] = 'G';
        ComplementBase['G'] = 'C';
        ComplementBase['T'] = 'A';
        return ComplementBase;
    }

    public:
    DNASeq() : seq(NULL), seq_len(0) {};
    DNASeq(const char* s) : seq(NULL), seq_len(0) {
        setSeq(s);
    }
    inline const char* getSeq() const {
        return seq;
    }
    inline void setSeq(const char* s) {
        if(s!=NULL) {
            seq = s;
            seq_len = strlen(seq);
        }
    }
    inline int size() const {
        return seq_len;
    }
    inline bool operator==(DNASeq const& a) const {
        return (a.seq != NULL) && (seq != NULL ) && (strcmp(seq, a.seq)) && (seq_len == a.seq_len);
    }
};

class Kmer {
    friend class Kmer_cmp;
    friend class Kmer_less;    
    friend class Kmer_hash;

    std::string kmer;
    unsigned int hash;

    void updateHash(){
        const int hash_basegroup = sizeof(unsigned int)*8/2;
        hash=0;
        for(size_t i=0; i<kmer.size(); i++) {
            switch(kmer[i]) {
                case 'A':
                    hash+=0<<((i%hash_basegroup)*2);
                    break;
                case 'C':
                    hash+=1<<((i%hash_basegroup)*2);
                    break;
                case 'G':
                    hash+=2<<((i%hash_basegroup)*2);
                    break;
                case 'T':
                    hash+=3<<((i%hash_basegroup)*2);
                    break;
                case 'Z':
                case 'N':
                    break;
                default:
                    throw std::runtime_error("Kmer::updateHash(): Kmer contains non-base character");
            }
        }
    }

    public:
    // Create Kmer from DNASeq and pos, kmer_len
    Kmer(std::tr1::shared_ptr<const DNASeq> seq, int pos, int len) {
        kmer.resize(len);
        for(int i=0;i<len; i++)
            kmer[i] = seq->seq[pos+i];
        updateHash();
    };

    // Create Kmer directly from string
    Kmer(std::string const& seq) {
        kmer = seq;
        updateHash();
    };

    inline std::string getKmer() const {
        return kmer;
    }
};

struct Kmer_hash {
    inline unsigned int operator()(Kmer const& kmer) const {
        return kmer.hash;
    }
};

struct Kmer_cmp {
    inline bool operator()(Kmer const& kmer1, Kmer const& kmer2) const {
        if(kmer1.hash!=kmer2.hash) return false;
        return (kmer1.kmer==kmer2.kmer);
    }
};

struct Kmer_less {
    inline bool operator()(Kmer const& kmer1, Kmer const& kmer2) const {
        if(kmer1.hash<kmer2.hash) return true;
        if(kmer1.hash>kmer2.hash) return false;
        return (kmer1.kmer<kmer2.kmer);
    }    
};

inline int baseToInt(char base) {
    switch(base) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'N':
        default:
            // Anything else will be regarded as 'N'
            return 4;
    }
    return 0;
};

inline char intToBase(int i) {
    switch(i) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            throw std::runtime_error("intToBase: argument must be in [0,3]");
    }
    return 0;
};

inline int intToBase(int i1, int i2) {
    if(i1==i2) return intToBase(i1);
    if(i1>i2) std::swap(i1, i2);
    if(i1==0) {
        if(i2==1) return 'M';
        if(i2==2) return 'R';
        if(i2==3) return 'W';
    } else if(i1==1) {
        if(i2==2) return 'S';
        if(i2==3) return 'Y';
    } else if(i1==2) {
        if(i2==3) return 'K';
    }
    return 0;
};

inline bool isAnyBase(int b){
    return (b>=4);
};

inline bool isAnyBase(char b){
    return (baseToInt(b)>=4);
};

#endif
