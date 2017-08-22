#define u32 unsigned int

// Number of input words.  Top 65 bits are overwritten with padding
// (512 is bits; 32 is bits-per-word.)
#define INPUT_SIZE	512 / 32  

// Number of words in a SHA-1 hash
#define OUTPUT_SIZE 160 / 32

// Number of words in a SHA-256 hash
#define OUTPUT_SIZE_256 256 / 32


struct shaInput {
	u32 msg[16];
};

struct shaOutput {
	u32 h[5];
};

struct shaOutput256 {
	u32 h[8];
};

u32 leftRotate(u32 val, int amount) {
	u32 r = (val << amount) | (val >> (32 - amount));
	return r;
}

u32 rightRotate(u32 val, int amount) {
	u32 r = (val >> amount) | (val << (32 - amount));
	return r;
}

u32 rightShift(u32 val, int amount) {
	u32 r = (val >> amount);
	return r;
}

u32 bit(u32 x, int i) {
	return (x & (1 << i)) >> i;
}

#define NUM_ROUNDS 80

// ---------------------- - variable - length SHA1

// handwritten multiplicative cost: 512 + (80-16)*32*3 + 80*(64+33) = 14416 
void sha1compress(u32 input[16], u32 hash[5])
{
	u32 W[80];  // Expanded message
	u32 i;
	u32 tmp;

	for (i = 0; i < 16; i += 1) {
		W[i] = input[i];
	}
	for (i = 16; i < NUM_ROUNDS; i += 1) {
		tmp = W[i - 3] ^ W[i - 8] ^ W[i - 14] ^ W[i - 16];
		W[i] = leftRotate(tmp, 1);
	}

	u32 a = hash[0];
	u32 b = hash[1];
	u32 c = hash[2];
	u32 d = hash[3];
	u32 e = hash[4];

	for (i = 0; i < NUM_ROUNDS; i += 1) {
		u32 f;
		u32 k;
		if (i < 20) {
#if PARAM==1
			// print(a);print(b);print(c);print(d);
			f = (b & c) | ((~b) & d); //CF costs only 32*x as  b*(c-d) + d 
#else
			f = 0;
			int j;
			for (j = 0; j < 32; j++)
				f += ((bit(b, j)*(bit(c, j) - bit(d, j)) + bit(d, j)) * (1 << j));
#endif
			k = 0x5A827999;
		}
		else if (i < 40) {
			f = b ^ c ^ d; // cost 32*2x
			//			f = 0;
			//			int j;
			//			for (j = 0; j < 32; j++) f += (bit(bit(b, j) + bit(c, j) + bit(d, j), 0) * (1 << j));
			k = 0x6ED9EBA1;
		}
		else if (i < 60) {
#if 0 
			f = (b & c) | (b & d) | (c & d); // cost 32*2x as we can use a binary decomposition: b+c+d = 2f + g
#else 
			f = 0;
			int j;
			for (j = 0; j < 32; j++) f += (bit(bit(b, j) + bit(c, j) + bit(d, j), 1) * (1 << j));
#endif
			k = 0x8F1BBCDC;
		}
		else if (i < 80) {
			f = b ^ c ^ d; // cost 32*2x
			//			f = 0;
			//			int j;
			//			for (j = 0; j < 32; j++) f += (bit(bit(b, j) + bit(c, j) + bit(d, j), 0) * (1 << j));
			k = 0xCA62C1D6;
		}
		tmp = leftRotate(a, 5) + f + e + k + W[i]; // free, but left with a 34-bit integer split in 34x 

		e = d;
		d = c;
		c = leftRotate(b, 30);
		b = a;
		a = tmp;
	}
	hash[0] += a;
	hash[1] += b;
	hash[2] += c;
	hash[3] += d;
	hash[4] += e;
}

// block-aligned 10* padding and length for both SHA1 and SHA256. 
void pad(u32 n, u32 msg[] /*input*/, u32 block[] /*output*/) {
	u32 bitlength = n * 8;
	u32 b = 1 + (bitlength + 8 + 63) / 512; // number of blocks
	u32 m = b * (512 / 32);                 // number of raw input words 

	// padding the message into b 512-bit blocks.
	u32 i;
	for (i = 0; 4 * i < n; i++) block[i] = msg[i];
	if (n % 4 == 0)
	{
		block[i++] = 0x80000000;;
	}
	else
	{
		block[i - 1] |= (0x80 << (24 - 8 * (n % 4)));
	}
	for (; i < m; i++)
	{
		block[i] = 0;
	}
	block[m - 1] |= ((bitlength / 256) << 8) + (bitlength % 256); // 64-bit big-endian length, < 2^32 of course
#if 0
	FILE* file;

	file = fopen("sha-input-padded-sha1n.bin", "wb");
	for (i = 0; i < m; i++)
	{
		unsigned char* charArr = (unsigned char*)&block[i];
		unsigned int reverse = charArr[3] << 0 | charArr[2] << 8
			| charArr[1] << 16 | charArr[0] << 24;
		fwrite(&reverse, sizeof(unsigned int), 1, file);
	}
	fclose(file);
#endif
}

void sha1n_lib(int* np /* number of bytes */, u32 * msg, struct shaOutput* digest) {
	u32 *hash = digest->h;
	int n = *np;

	u32 bitlength = n * 8;
	u32 b = 1 + (bitlength + 8 + 63) / 512; // number of blocks
	u32 m = b * (512 / 32);                 // number of raw input words 

	// padding the message into b 512-bit blocks.

#ifdef MQAP
	u32 input[m];
#else
	//u32 input[(LENGTH + 3) / 4 + 16]; //MK: only in QAP we know np at runtime.
	u32* input;
	input = (u32*) malloc(m*4);
#endif 
	u32 i;
	pad(n, msg, input);

	hash[0] = 0x67452301;
	hash[1] = 0xEFCDAB89;
	hash[2] = 0x98BADCFE;
	hash[3] = 0x10325476;
	hash[4] = 0xC3D2E1F0;

	for (i = 0; i<b; i++) sha1compress(input + (i*(512 / 32)), hash);
}

#define NUM_ROUNDS_256 64

void chunk256(u32 input[16], u32 hash[8])
{
	u32 W[64];  // Expanded message
	u32 i;
	u32 s0;
	u32 s1;

	/*
	u32 K[64] = { 0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
	0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
	0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
	0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
	0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
	0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
	0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
	0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2 };
	*/

	u32 K[64]; // to avoid global constants...
	i = 0;
	K[i++] = 0x428a2f98;
	K[i++] = 0x71374491;
	K[i++] = 0xb5c0fbcf;
	K[i++] = 0xe9b5dba5;
	K[i++] = 0x3956c25b;
	K[i++] = 0x59f111f1;
	K[i++] = 0x923f82a4;
	K[i++] = 0xab1c5ed5;

	K[i++] = 0xd807aa98;
	K[i++] = 0x12835b01;
	K[i++] = 0x243185be;
	K[i++] = 0x550c7dc3;
	K[i++] = 0x72be5d74;
	K[i++] = 0x80deb1fe;
	K[i++] = 0x9bdc06a7;
	K[i++] = 0xc19bf174;

	K[i++] = 0xe49b69c1;
	K[i++] = 0xefbe4786;
	K[i++] = 0x0fc19dc6;
	K[i++] = 0x240ca1cc;
	K[i++] = 0x2de92c6f;
	K[i++] = 0x4a7484aa;
	K[i++] = 0x5cb0a9dc;
	K[i++] = 0x76f988da;

	K[i++] = 0x983e5152;
	K[i++] = 0xa831c66d;
	K[i++] = 0xb00327c8;
	K[i++] = 0xbf597fc7;
	K[i++] = 0xc6e00bf3;
	K[i++] = 0xd5a79147;
	K[i++] = 0x06ca6351;
	K[i++] = 0x14292967;

	K[i++] = 0x27b70a85;
	K[i++] = 0x2e1b2138;
	K[i++] = 0x4d2c6dfc;
	K[i++] = 0x53380d13;
	K[i++] = 0x650a7354;
	K[i++] = 0x766a0abb;
	K[i++] = 0x81c2c92e;
	K[i++] = 0x92722c85;

	K[i++] = 0xa2bfe8a1;
	K[i++] = 0xa81a664b;
	K[i++] = 0xc24b8b70;
	K[i++] = 0xc76c51a3;
	K[i++] = 0xd192e819;
	K[i++] = 0xd6990624;
	K[i++] = 0xf40e3585;
	K[i++] = 0x106aa070;

	K[i++] = 0x19a4c116;
	K[i++] = 0x1e376c08;
	K[i++] = 0x2748774c;
	K[i++] = 0x34b0bcb5;
	K[i++] = 0x391c0cb3;
	K[i++] = 0x4ed8aa4a;
	K[i++] = 0x5b9cca4f;
	K[i++] = 0x682e6ff3;

	K[i++] = 0x748f82ee;
	K[i++] = 0x78a5636f;
	K[i++] = 0x84c87814;
	K[i++] = 0x8cc70208;
	K[i++] = 0x90befffa;
	K[i++] = 0xa4506ceb;
	K[i++] = 0xbef9a3f7;
	K[i++] = 0xc67178f2;

	for (i = 0; i < 16; i += 1) {
		W[i] = input[i];
	}
	for (i = 16; i < NUM_ROUNDS_256; i += 1) {
		s0 = rightRotate(W[i - 15], 7) ^ rightRotate(W[i - 15], 18) ^ rightShift(W[i - 15], 3);
		s1 = rightRotate(W[i - 2], 17) ^ rightRotate(W[i - 2], 19) ^ rightShift(W[i - 2], 10);
		W[i] = W[i - 16] + s0 + W[i - 7] + s1;
	}
	u32 a = hash[0];
	u32 b = hash[1];
	u32 c = hash[2];
	u32 d = hash[3];
	u32 e = hash[4];
	u32 f = hash[5];
	u32 g = hash[6];
	u32 h = hash[7];

	u32 S1;
	u32 ch;
	u32 temp1;
	u32 S0;
	u32 maj;
	u32 temp2;

	for (i = 0; i < NUM_ROUNDS_256; i += 1) {
		S1 = rightRotate(e, 6) ^ rightRotate(e, 11) ^ rightRotate(e, 25);
		ch = (e & f) ^ ((~e) & g);
		temp1 = h + S1 + ch + K[i] + W[i];
		S0 = rightRotate(a, 2) ^ rightRotate(a, 13) ^ rightRotate(a, 22);
		maj = (a & b) ^ (a & c) ^ (b & c);
		temp2 = S0 + maj;

		h = g;
		g = f;
		f = e;
		e = d + temp1;
		d = c;
		c = b;
		b = a;
		a = temp1 + temp2;
	}
	hash[0] += a;
	hash[1] += b;
	hash[2] += c;
	hash[3] += d;
	hash[4] += e;
	hash[5] += f;
	hash[6] += g;
	hash[7] += h;
}

int sha256n_lib(int* np /* number of bytes */, u32 *msg, struct shaOutput256* digest) {
	u32 *hash = digest->h;
	int n = *np;

	u32 bitlength = n * 8;
	u32 b = 1 + (bitlength + 8 + 63) / 512; // number of blocks
	u32 m = b * (512 / 32);                 // number of raw input words 

	// padding the message into b 512-bit blocks.
    #ifdef MQAP
	u32 input[m];
    #else
	u32* input;
	input = (u32*)malloc(m * 4);
	//u32 input[(LENGTH + 3) / 4 + 16];
    #endif 

	u32 i;

	pad(n, msg, input);

	hash[0] = 0x6a09e667;
	hash[1] = 0xbb67ae85;
	hash[2] = 0x3c6ef372;
	hash[3] = 0xa54ff53a;
	hash[4] = 0x510e527f;
	hash[5] = 0x9b05688c;
	hash[6] = 0x1f83d9ab;
	hash[7] = 0x5be0cd19;

	for (i = 0; i<b; i++) chunk256(input + (i*(512 / 32)), hash);
	return 0;
}
