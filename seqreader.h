
/**
 * @file seqreader.h
 *
 * @brief sequence reader and alignment writers
 * 
 * @detail
 * Sequence readers and alignment writers are provided to abstract
 * read/write operation to sequences. This abstraction enables the
 * alignment implementations in this library to support different 
 * input/output formats with the same source code.
 * 
 * usage:
 * 
 * DECLARE_SEQ(name) declares a pointer to the sequence array and
 * read buffer.
 * CLEAR_SEQ(name, ptr) inirializes pointer to sequence.
 *
 * FETCH(name, pos) fetches (reads) a base at position pos. the fetched
 * base will be used in the following DECODE(name) or COMPARE(name1, name2)
 *
 *
 * To specify the input/output formats, you need to give proper 
 * defines to the compiler with -D option, such as:
 * 
 * $ gcc -DSEQ=seq8 -DALN=str8 ...
 *
 * the full list of the required defines and options are:
 *
 * a) Input sequence format. specified by SEQ definition
 *
 * ascii:
 * sequences represented in ascii string. a string like "ATAA".
 * only capital cases are allowed.
 *
 *
 * seq4:
 * sequences are stored in an array of 8-bit unsigned. each element
 * consists a bit-field of bases. at least one of the four field
 * must be set.
 *
 * the configuration is following:
 * | bit 7  ...........  4 | 3 | 2 | 1 | 0 |
 * | unused (must be zero) | T | C | G | A |
 *
 *
 * seq2:
 * sequences are stored in an array of 8-bit unsigned. each element
 * holds bases in the lower 2-bit. this is a native format of this 
 * library.
 *
 * the corresopondence are following:
 * A: 0x00
 * C: 0x01
 * G: 0x02
 * T: 0x03
 *
 *
 * seq4p8:
 * sequences are stored in an array of 8-bit unsigned. each element
 * hold two bases in it. the lower 4-bit represents earlier base,
 * the upper 4-bit represents following base. each 4-bits consists
 * a bit-field of {A, C, G, T}, the same as 'seq4'.
 *
 * seq2p8:
 *
 * seq1p64:
 * seqences are stored in an array of 64-bit unsigned. each 64-bit
 * contains a series of upper/lower 1-bit of 'seq2' 2-bit format.
 * 
 *
 *
 * b) output alignment format
 * cigar:
 * returns a ASCII cigar string. "M" stands for matches, "X" for mismatches,
 * "I" for insertion to sequence a, and "D" for deletions from sequence a.
 * a simple examples are:
 *   a = "AAA" and b = "AAA" yields "3M"
 *   a = "AAAT" and b = "AAT" yields "2M1D1M" ("1D" means a deletion of "A" from sequence a)
 *
 * ascii:
 * returns an string of {'M', 'X', 'I', 'D'}. each character represents
 * a match of bases in two sequences, a mismatch of bases, an insertion of
 * a base to sequence a, and a deletion of a base from sequence b.
 * a simple example are:
 *   a = "AAA" and b = "AAA" yields "MMM"
 *   a = "AAAT" and b = "AAT" yields "MMDM"
 *
 * see also: wscript in each directory, alnbld.py in util, table.th in util.
 */

#ifndef _SEQREADER_H_INCLUDED
#define _SEQREADER_H_INCLUDED

#include <stdio.h>
#include <string.h>						/* for memmove */
#include "sea.h"				/* for definitions of sea_int_t */

/**
 * labels
 */
#define ascii 				( 0 )
#define seq4 				( 1 )
#define seq2 				( 2 )
#define seq4p8 				( 3 )
#define seq2p8				( 4 )
#define seq1p64				( 5 )
#define DIM_SEQ		 		( 6 ) 		/* equals to seq1p64+1 */

#define ascii 				( 0 )
#define cigar 				( 1 )
#define DIM_ALN				( 2 ) 		/* equals to cig16+1 */


/**
 * check SEQ and ALN is defined, and the value is in 0..DIM_SEQ or 0..DIM_ALN
 */
#define SEQ 				seq2
#define ALN 				cigar

#ifndef SEQ
 	#error "SEQ must be defined"
#endif

#if SEQ < 0 || SEQ >= DIM_SEQ
 	#error "SEQ must satisfy 0 < SEQ < "##DIM_SEQ
#endif

#ifndef ALN
 	#error "ALN must be defined"
#endif

#if ALN < 0 || ALN >= DIM_ALN
 	#error "ALN must satisfy 0 < ALN < "##DIM_ALN
#endif

/**
 * input sequence reader macros.
 *
 * declearation and initialization of a pointer.
 * _##name##seqptr: a pointer to the head of the input sequence.
 * _##name##seqbuf: read buffer
 */
#if SEQ == seq1p64
	#define DECLARE_SEQ(name)		unsigned long const *_##name##seqptr; \
 									unsigned long _##name##seqbufl, _##name##seqbufh; \
 									unsigned long _##name##spos;
	#define CLEAR_SEQ(name,base,pos,len) \
									_##name##seqptr = (unsigned long const *)(base); \
									_##name##spos = pos;

#else
	#define DECLARE_SEQ(name)		char const *_##name##seqptr; \
 									char _##name##seqbuf; \
 									unsigned long _##name##spos;
	#define CLEAR_SEQ(name,base,pos,len) \
									_##name##seqptr = (char const *)(base); \
									_##name##spos = pos;

#endif

/**
 * definitions of a random access reader.
 */
#if SEQ == ascii
 	/*
 	 * ASCII character to 'seq2' conversion.
 	 * the bit 1 and the bit 2 of {'A', 'C', 'G', 'T'} or {'a', 'c', 'g', 't'}
 	 * encodes four different numbers, the conversion from ASCII to 'seq2'
 	 * result in just bit operations between this two bits.
 	 * 
 	 * when an ASCII character a is given, an actual operation forms:
 	 * ((a>>1) ^ (a>>2)) & 0x03
 	 */

 	/**
 	 * @macro FETCH
 	 * @brief issue a random access to sequence.
 	 * a < x && x < b x - a < b - a
 	 */
 	#define FETCH(name, pos)		{ _##name##seqbuf = _##name##seqptr[_##name##spos + pos]; }
 	
 	/**
 	 * @macro DECODE
 	 * @brief decode a fetched ascii character (a base) to seq2 number.
 	 */
 	#define DECODE(name)			( _##name##seqbuf )
// 	#define DECODE(name) 			( ((_##name##seqbuf>>1) ^ (_##name##seqbuf>>2)) & 0x03 )
 	
 	/**
 	 * @macro DECODE_LOW
 	 * @brief extract the lower bit of seq2 number.
 	 */
 	#define DECODE_LOW(name)	 	( ((_##name##seqbuf>>1) ^ (_##name##seqbuf>>2)) & 0x01 )
 	
 	/**
 	 * @macro DECODE_HIGH
 	 * @brief extract the higher bit of seq2 number.
 	 */
 	#define DECODE_HIGH(name) 		( ((_##name##seqbuf>>2) ^ (_##name##seqbuf>>3)) & 0x01 )
 	
 	/**
 	 * @macro COMPARE
 	 * @brief compare two fetched characters.
 	 */
 	#define COMPARE(name1, name2)	( _##name1##seqbuf == _##name2##seqbuf )

#elif SEQ == seq4
 	/*
 	 * seq4
 	 */
 	#define FETCH(name, pos)		{ _##name##seqbuf = _##name##seqptr[_##name##spos + pos]; }
 	#define DECODE(name)			( _##name##seqbuf )
// 	#define DECODE(name)			( ((_##name##seqbuf>>1) - ((_##name##seqbuf>>3) & 0x01)) & 0x03 )
 	#define DECODE_LOW(name) 		( ((_##name##seqbuf>>1) | ((_##name##seqbuf>>3) & 0x01)) )
 	#define DECODE_HIGH(name) 		( ((_##name##seqbuf>>2) | ((_##name##seqbuf>>3) & 0x01)) )
 	#define COMPARE(name1, name2)	( _##name1##seqbuf == _##name2##seqbuf )

#elif SEQ == seq2
 	/*
 	 * seq2 random accessor
 	 */
 	#define FETCH(name, pos)		{ _##name##seqbuf = _##name##seqptr[_##name##spos + pos]; }
 	#define DECODE(name)			( _##name##seqbuf )
 	#define DECODE_LOW(name) 		( _##name##seqbuf & 0x01 )
 	#define DECODE_HIGH(name) 		( (_##name##seqbuf>>1) & 0x01 )
 	#define COMPARE(name1, name2)	( _##name1##seqbuf == _##name2##seqbuf )

#elif SEQ == seq4p8
 	/*
 	 * 8-bit packed seq4
 	 */
 	#define FETCH(name, pos)		{ _##name##seqbuf = _##name##seqptr[(_##name##spos + pos)/2]>>(((_##name##spos + pos)&0x01)<<2); }
 	#define DECODE(name)			( _##name##seqbuf & 0x0f )
// 	#define DECODE(name)			( ((_##name##seqbuf>>1) - ((_##name##seqbuf>>3) & 0x01)) & 0x03 )
 	#define DECODE_LOW(name) 		( ((_##name##seqbuf>>1) | ((_##name##seqbuf>>3) & 0x01)) )
 	#define DECODE_HIGH(name) 		( ((_##name##seqbuf>>2) | ((_##name##seqbuf>>3) & 0x01)) )
 	#define COMPARE(name1, name2)	( _##name1##seqbuf == _##name2##seqbuf )

#elif SEQ == seq2p8
 	/*
 	 * 8-bit packed seq2
 	 */
 	#define FETCH(name, pos)		{ _##name##seqbuf = _##name##seqptr[(_##name##spos + pos)/4]>>(((_##name##spos + pos)&0x02)<<1); }
 	#define DECODE(name)			( _##name##seqbuf & 0x03 )
 	#define DECODE_LOW(name) 		( _##name##seqbuf & 0x01 )
 	#define DECODE_HIGH(name) 		( (_##name##seqbuf>>1) & 0x01 )
 	#define COMPARE(name1, name2)	( _##name1##seqbuf == _##name2##seqbuf )

#elif SEQ == seq1p64
 	/*
 	 * seq1p64
 	 */
 	#define FETCH(name, pos)		{ _##name##seqbufl = _##name##seqptr[((_##name##spos + pos)/32) & ~0x01]>>((_##name##spos + pos)&0x3f); \
 									  _##name##seqbufh = _##name##seqptr[((_##name##spos + pos)/32) | 0x01]>>((_##name##spos + pos)&0x3f); }
 	#define DECODE(name)			( _##name##seqbufl & 0x01 | ((_##name##seqbufh & 0x01)<<1) )
 	#define DECODE_LOW(name)		( _##name##seqbufl & 0x01 )
 	#define DECODE_HIGH(name)		( _##name##seqbufh & 0x01 )
 	#define COMPARE(name1, name2)	( _##name1##seqbufl == _##name2##seqbufl && _##name1##seqbufh == _##name2##seqbufh )

#else
 	#error "invalid value of SEQ"
#endif

#if 0
/**
 * alignment writer macros
 *
 * DECLARE_ALN:
 * declares a pointer to the head of the buffer (_##name##alnbase),
 * a pointer to a current tail of the buffer (_##name##alnptr),
 * and a pointer to the end of the buffer (_##name#sent) (work as a sentinel).
 *
 * CLEAR_ALN: initialize _##name##alnbase, _##name##alnptr, and _##name##sent.
 *
 * PUSH: push a character to the tail of the buffer and increment the pointer.
 *
 * SWAP_CHAR: swap two characters with temporary register (tmp)
 *
 * REVERSE: reverse a content of the buffer. (must be called once when the
 *          traceback is finished.)
 *
 * LENGTH: returns the length of the alignment string. this function MUST NOT
 *         be called once the REVERSE function is executed.
 */
#if ALN == ascii
	/*
	 * ascii string
	 */
	#define MATCH_CHAR 				'M'
	#define MISMATCH_CHAR			'X'
	#define INSERTION_CHAR 			'I'
	#define DELETION_CHAR			'D'
/*
 	#define DECLARE_ALN(name) 		char *_##name##alnptr, *_##name##alnbase, *_##name##sent;
	#define CLEAR_ALN(name, ptr, len) { \
		_##name##alnbase = ptr; \
		_##name##alnptr = (char *)ptr + len; \
		_##name##sent = (char *)ptr + len; \
	}
	#define PUSH(name, c) { \
		_##name##alnptr -= (_##name##alnptr > _##name##alnbase); \
		*_##name##alnptr = (c); \
	}
	#define REVERSE(name, tmp) { \
		while(_##name##alnptr < _##name##sent) { \
			*_##name##alnbase++ = *_##name##alnptr++; \
		} \
		*_##name##alnbase = '\0'; \
	}
	#define LENGTH(name) 			( _##name##sent - _##name##alnptr )
*/

	#define DECLARE_ALN(name) 		char *_##name##alnptr, *_##name##alnbase;
	
	#define CLEAR_ALN(name, ptr, len) { \
		_##name##alnbase = ptr; \
		_##name##alnptr = ptr; \
	}

	#define PUSH(name, c) { \
		*_##name##alnptr++ = (c); \
	}

	#define SWAP_CHAR(a, b, tmp) { \
		(tmp) = (a); (a) = (b); (b) = (char)(tmp); \
	}

	#define REVERSE(name) { \
		char tmp; \
		*_##name##alnptr-- = '\0'; \
		while(_##name##alnbase < _##name##alnptr) { \
			SWAP_CHAR(*_##name##alnbase, *_##name##alnptr, tmp); \
			_##name##alnbase++; _##name##alnptr--; \
		} \
	}

	#define LENGTH(name) 			( _##name##alnptr - _##name##alnbase )

#elif ALN == cigar
 	/*
 	 * cigar string
 	 */
	#define MATCH_CHAR 				'M'
	#define MISMATCH_CHAR			'X'
	#define INSERTION_CHAR 			'I'
	#define DELETION_CHAR			'D'
 	#define DECLARE_ALN(name) 		char *_##name##alnptr, *_##name##alnbase, *_##name##tmpptr; \
 	 								char _##name##c, _##name##buf[16]; \
 	 								unsigned long _##name##cnt, _##name##l;
	
	#define CLEAR_ALN(name, ptr, len) { \
		_##name##alnbase = ptr; \
		_##name##alnptr = ptr + len; \
		_##name##c = 'R'; \
		_##name##cnt = 0; \
		_##name##l = 0; \
	}

	#define PUSH(name, ch) { \
		_##name##l++; \
		if((ch) == _##name##c) { \
			_##name##cnt++; \
		} else { \
			_##name##tmpptr = _##name##buf; \
			sprintf(_##name##tmpptr, "%lu", ++_##name##cnt); \
			while(*_##name##tmpptr != '\0') { *_##name##tmpptr++; } \
			*_##name##alnptr++ = _##name##c; \
			while(_##name##tmpptr != _##name##buf) { \
				*_##name##alnptr++ = *--_##name##tmpptr; \
			} \
			_##name##c = (ch); \
			_##name##cnt = 0; \
		} \
	}

	#define SWAP_CHAR(a, b, tmp) { \
		(tmp) = (a); (a) = (b); (b) = (char)(tmp); \
	}
	#define REVERSE(name) { \
		char tmp; \
		_##name##tmpptr = _##name##buf; \
		sprintf(_##name##tmpptr, "%lu", ++_##name##cnt); \
		while(*_##name##tmpptr != '\0') { *_##name##tmpptr++; } \
		*_##name##alnptr++ = _##name##c; \
		while(_##name##tmpptr != _##name##buf) { \
			*_##name##alnptr++ = *--_##name##tmpptr; \
		} \
		_##name##cnt = 0; \
		_##name##tmpptr = --_##name##alnptr; \
		while(_##name##alnbase < _##name##alnptr) { \
			SWAP_CHAR(*_##name##alnbase, *_##name##alnptr, tmp); _##name##alnbase++; _##name##alnptr--; \
		} \
		*(_##name##tmpptr-1) = '\0'; \
	}

	#define LENGTH(name) 			( _##name##l + _##name##cnt )
#else
 	#error "invalid value of ALN"
#endif
#endif

#endif /* #ifndef _SEQREADER_H_INCLUDED */

/*
 * end of seqreader.h
 */
