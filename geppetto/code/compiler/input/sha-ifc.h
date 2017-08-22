#if PARAM==0
#elif PARAM==1
#elif PARAM==2
#else
#error unknown PARAM
#endif



struct bank_Input {
//	u32 msg[INPUT_SIZE];
	struct shaInput i;
};

struct bank_Output {
//	u32 h[OUTPUT_SIZE];
	struct shaOutput o;
};

struct bank_Output256 {
//	u32 h[OUTPUT_SIZE_256];
	struct shaOutput256 o;
};

//void outsource(struct bank_Input*, struct bank_Output*);

//--- generalized tests, for a parametric length (define LENGTH in chars)


struct bank_Msg {
	struct shaMsg m;
};
