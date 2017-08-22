#ifdef WHOLE
#define NUM_BANKS 2
#include "geppetto.h"
#endif

#include "pinocchio.h"
#ifndef _WIN32
#define MQAP 1
#endif


#ifndef MQAP
#include <stdio.h>
#include <string.h>
#include "pinocchio.c"
#endif

// Assumes the message is exactly 416 bits long, 
// so the value we hash is: msg || 1000...00 || len,
// where len = 416 and there are 31 zeros
// This saves us from dealing with padding for now

// Based in part on pseudocode from http://en.wikipedia.org/wiki/SHA-1
// the spec at http://tools.ietf.org/html/rfc3174,
// and the FIPS 180-2 spec: 
// http://csrc.nist.gov/publications/fips/fips180-2/fips180-2.pdf


// PARAM==1 : normal 
// PARAM==2 : with binary operations optimized for QAPs (not so diffferent)

#include "sha-library.c"

#ifndef LENGTH
#define LENGTH 3
#endif

struct shaMsg {
	u32 msg[(LENGTH + 3) / 4];
};

#include "sha-ifc.h"

// handwritten multiplicative cost: 512 + (80-16)*32*3 + 80*(64+33) = 14416 
void outsource(struct bank_Input *input0, struct bank_Output *output)
{
#ifdef WRITE_FILE
	FILE* file;
#endif 
	u32* hash;
	u32 W[80];  // Expanded message
	u32 i;
	u32 tmp;
	u32 a, b, c, d, e;


	// making a local copy, so that input0 is input-only
	struct shaInput input[1];
	for (i = 0; i < INPUT_SIZE; i++) {
		input->msg[i] = input0->i.msg[i];
	}

	// a hack, ensuring that the last 3 words of the input are read 
	// (they are equal to 0 in the tests)
	input->msg[0] += input->msg[INPUT_SIZE - 1];
	input->msg[0] += input->msg[INPUT_SIZE - 2];
	input->msg[0] += input->msg[INPUT_SIZE - 3];

	input->msg[INPUT_SIZE - 1] = 0x000001A0;
	input->msg[INPUT_SIZE - 2] = 0;          // (length is 64 bits)
	input->msg[INPUT_SIZE - 3] = 0x80000000; // Fix the padding

#ifdef WRITE_FILE
	// Write out padded input 
	file = fopen("sha-input-padded.bin", "wb");
	for (i = 0; i < INPUT_SIZE; i += 1) {
		unsigned char* charArr = (unsigned char*)&input->msg[i];
		int reverse = charArr[3] << 0 | charArr[2] << 8 |
			charArr[1] << 16 | charArr[0] << 24;
		fwrite(&reverse, sizeof(int), 1, file);
	}
	//fwrite(input->msg, sizeof(char), (INPUT_SIZE)*sizeof(u32), file); 
	//fwrite(input->msg, sizeof(u32), INPUT_SIZE, file); //MK commented out
	fclose(file);
#endif
	for (i = 0; i < 16; i += 1) {
		W[i] = input->msg[i];
		//print(W[i]);
	}
	// Entry cost: decompose 512-bit input
	// Expand the message further
	for (i = 16; i < NUM_ROUNDS; i += 1) {
		u32 tmp = W[i - 3] ^ W[i - 8] ^ W[i - 14] ^ W[i - 16]; // cost: (80-16)*32*3 x except for constant propagations
		W[i] = leftRotate(tmp, 1);
		//print(W[i]);
	}

	hash = output->o.h;
	// Initialize the state variables
	hash[0] = 0x67452301;
	hash[1] = 0xEFCDAB89;
	hash[2] = 0x98BADCFE;
	hash[3] = 0x10325476;
	hash[4] = 0xC3D2E1F0;

	a = hash[0];
	b = hash[1];
	c = hash[2];
	d = hash[3];
	e = hash[4];

	//#define TICKER
#ifdef TICKER
#define TICK log0 = log; log = nRoot(); print(log - log0)
	int log0, log = 0;
#else 
#define TICK {}
#endif

	// hoping for 97+ bit multiplications per iteration. 
	for (i = 0; i < NUM_ROUNDS; i += 1) {
		u32 f;
		u32 k;
		//		TICK;

		if (i < 20) {
#if PARAM==1
			// print(a);print(b);print(c);print(d);
			f = (b & c) | ((~b) & d); //CF costs only 32* as  b*(c-d) + d 
#else
			f = 0;
			int j;
			for (j = 0; j < 32; j++)
				f += bound(bit(b, j)*(bit(c, j) - bit(d, j)) + bit(d, j), 0, 1) * (1 << j);

			//print(f);
#endif
			k = 0x5A827999;
		}
		else if (i < 40) {
#if 1
			f = b ^ c ^ d; // cost 32*2x
#else
			f = 0;
			int j;
			for (j = 0; j < 32; j++) f += bit(bit(b, j) + bit(c, j) + bit(d, j), 0) * (1 << j);
#endif
			k = 0x6ED9EBA1;
		}
		else if (i < 60) {
#if PARAM==1 
			f = (b & c) | (b & d) | (c & d); // cost 32*2x as we can use a binary decomposition: b+c+d = 2f + g
#else 
			f = 0;
			int j;
			for (j = 0; j < 32; j++) f += bit(bit(b, j) + bit(c, j) + bit(d, j), 1) * (1 << j);
#endif
			k = 0x8F1BBCDC;
		}
		else if (i < 80) {
#if 1
			f = b ^ c ^ d; // cost 32*2x
#else
			f = 0;
			int j;
			for (j = 0; j < 32; j++) f += (bit(bit(b, j) + bit(c, j) + bit(d, j), 0) * (1 << j));
#endif
			k = 0xCA62C1D6;
		}

		//		TICK;
		tmp = leftRotate(a, 5) + f + e + k + W[i]; // free, but left with a 35-bit integer split. 

		//		TICK;
		e = d;
		d = c;
		c = leftRotate(b, 30);
		b = a;
		a = tmp;

#ifndef MQAP 
		// printf("round %d (f,c,a).\n", i);
#endif
		// the new values getting in.
		// print(f); print(c);	print(a);
	}

	hash[0] = hash[0] + a;
	hash[1] = hash[1] + b;
	hash[2] = hash[2] + c;
	hash[3] = hash[3] + d;
	hash[4] = hash[4] + e;
}


void sha1n(int* np /* number of bytes */, struct bank_Msg* message, struct bank_Output* digest) {
	sha1n_lib(np, message->m.msg, &(digest->o));
}

void sha256n(int* np /* number of bytes */, struct bank_Msg* message, struct bank_Output256* digest) {
	sha256n_lib(np, message->m.msg, &(digest->o));
}

// variant: how to hide n? 

#ifndef MQAP

void print_ints(int size, int ints[]) {
	int i;
	for (i = 0; i < size; i++) {
		printf("%08x ", ints[i]);
	}
//	printf("\n");
}

void bytes_to_ints(char chars[], u32 ints[]) {
	int i;
	int n = (int) strlen(chars);
	for (i = 0; i < n; i++){
		ints[i / 4] = 0;
	}
	for (i = 0; i < n; i++){
		ints[i / 4] = ints[i / 4] << 8;
		ints[i / 4] = ints[i / 4] + chars[i];
	}
	int j; //shift the last word correctly
	for (j = 0; j < (4-(n % 4)) % 4; j++) {
		ints[(n-1)/4] = ints[(n-1)/4] << 8;
	}
}

void quick_brown_fox_test() {
	char m[] = "The quick brown fox jumps over the lazy dog";
	printf("Wikipedia example: '%s', length %d \n", m, strlen(m));

	int length = 43;

	u32 msg[(43 + 3) / 4];

	bytes_to_ints(m, msg);
	printf("The input in ints is : \n");
	print_ints((length + 4) / 4, msg);
	printf("\n");

	struct shaOutput shaOut;

	sha1n_lib(&length, msg, &shaOut);
	printf("The output for quick fox is : ", m);
	print_ints(OUTPUT_SIZE, shaOut.h);
	printf("\nThe output should be        : 2fd4e1c6 7a2d28fc ed849ee1 bb76e739 1b93eb12\n");
}

void variable_size_test() {
	// test variable-size
	int length = 4;
	u32 msg[1]; //needed because bytes_to_ints adds a 0 to end the string.
	struct shaOutput shaOut;

	bytes_to_ints("abc ", msg);
	printf("\nThe input for 'abc ' in ints is: ");
	print_ints(1, msg);
	printf("\nThe input for 'abc ' should be : 61626320\n");

	
	sha1n_lib(&length, msg, &shaOut);
	printf("\nThe sha1   output for 'abc ' is: ");
	print_ints(OUTPUT_SIZE, shaOut.h);
	printf("\n");

	struct shaOutput256 shaOut256;
	int n;
	n = 3;
	sha256n_lib(&n, msg, &shaOut256);
	printf("The sha256 output for 'abc ' is: ");
	print_ints(OUTPUT_SIZE, shaOut256.h);
	printf("\n");
}

void iterative_outsource_test() {
	struct bank_Input  input;
	struct bank_Output output;
	int i, iter;

	// independent iterations; cut to 1 for debugging.
	for (iter = 4; iter < 5; iter++) {

		for (i = 0; i<INPUT_SIZE - 3; i++)
		{
			input.i.msg[i] = (i + iter) << (iter * 3);
		}
		for (i = INPUT_SIZE - 3; i<INPUT_SIZE; i++)
		{
			input.i.msg[i] = 0;
		}


		//original padding for outsource
		for (i = 1; i < 512 / 32; i++) input.i.msg[i] = 0;

#ifdef WRITE_FILE
		// Write out a file we can test with sha1sum
		char *sha_in_fn = "win-input.bin";
		FILE* file;

		file = fopen(sha_in_fn, "wb");
		//fwrite("abc", sizeof(char), 3, file);
		for (i = 0; i < INPUT_SIZE - 3; i++)
		{
			unsigned char* charArr = (unsigned char*)&input.msg[i];
			unsigned int reverse = charArr[3] << 0 | charArr[2] << 8
				| charArr[1] << 16 | charArr[0] << 24;
			//printf("%x -> %x, %x %x %x %x\n", input->msg[i], reverse, charArr[3], charArr[2], charArr[1], charArr[0]);
			fwrite(&reverse, sizeof(unsigned int), 1, file);
		}
		fclose(file);
#endif
		outsource(&input, &output);

		printf("The output for '");
		print_ints(INPUT_SIZE, input.i.msg);
		printf("' through outsource is: ");
		print_ints(OUTPUT_SIZE, output.o.h);
		printf("\n");
	}
}

int main() {
#ifdef WHOLE
	init();
#endif
	quick_brown_fox_test();
	variable_size_test();
	iterative_outsource_test();
	// getc(stdin);
	return 0;

}
#endif
