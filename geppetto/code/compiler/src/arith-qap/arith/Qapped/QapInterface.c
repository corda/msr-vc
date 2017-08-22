#include "QapInterface.h"

#include "QapFp4.h"
#include "QapFp6.h"
#include "QapFp12.h"
#include "Qappairing.h"
#include "Qapmiller.h"
#include "Qapfinalexpo.h"

void QapInterface_init(QapInterface* ifc)
{
	  QapFp_initialize_config();       
  	QapFp_init(ifc->m);
#ifndef TINY
  /* This is the 254-bit BN prime p = 3 (mod 4). */
  /* A point on the BN curve is composed of two values mod this prime */
	QapFp_set_u32(
	  ifc->m,
	  0x25236482, 0x40000001, 0xBA344D80, 0x00000008,
	  0x61210000, 0x00000013, 0xA7000000, 0x00000013);
#else
    /* This is the 53-bit BN prime p = 3 (mod 4). */
    /* A point on the BN curve is composed of two values mod this prime */
    QapFp_set_u32(ifc->m, 0, 0, 0, 0, 0, 0, 0x15C9B0, 0x796B71DB);
#endif // TINY
  QapFp_init(ifc->b);
	QapFp_wp_init (ifc->Lg);
	QapFp2_initialize_config (); 

	QapFp2_init (ifc->B);
	QapFp2_init (ifc->B3);
	QapFp2_init (ifc->l00);
	QapFp2_init (ifc->l01);
	QapFp2_init (ifc->l10);
	
	QapFp_set_ui (ifc->b, 2);
	
#ifndef TINY
  QapFp_set_u32(ifc->Lg->X,
    0x25236482, 0x40000001, 0xBA344D80, 0x00000008,
    0x61210000, 0x00000013, 0xA7000000, 0x00000012 );

  QapFp_set_u32(ifc->Lg->Y,
    0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000001 );

	QapFp_set_u32(ifc->B->a1,
		0x25236482, 0x40000001, 0xBA344D80, 0x00000008,
		0x61210000, 0x00000013, 0xA7000000, 0x00000012);
	QapFp_set_u32(ifc->B->a0,
		0x00000000, 0x00000000, 0x00000000, 0x00000000,
		0x00000000, 0x00000000, 0x00000000, 0x00000001);
#else
  QapFp_set_u32(ifc->Lg->X,
    0, 0, 0, 0,
    0, 0, 0x32134, 0x7FC84CEE);

  QapFp_set_u32(ifc->Lg->Y,
    0, 0, 0, 0,
    0, 0, 0xC5D3F, 0x93EEC09B);

  QapFp_set_u32(ifc->B->a1,
    0, 0, 0, 0,
    0, 0, 0x15C9B0, 0x796B71DA);
  QapFp_set_u32(ifc->B->a0,
    0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000001);
#endif //TINY
    QapFp_initialize_weierstrass (&ifc->b, 1);

	QapFp2_initialize_weierstrass (ifc->B);

	QapFp2_wp_init (ifc->Rg);

	QapFp4_initialize_config ();
	QapFp6_initialize_config ();
	QapFp12_initialize_config ();
	
	Qapbn_bls12_miller_initialize_config (ifc->B);
	Qapbn_bls12_finalexpo_initialize_config ();


#ifndef TINY
	// Yes, the sample code really does load these in reverse array
	// index order; I guess because that way everything is little-endian.
	// I LOOOVE little-endian!
	QapFp_set_u32(ifc->Rg->X->a1,
    0x0516AAF9, 0xBA737833, 0x310AA78C, 0x5982AA5B,
    0x1F4D746B, 0xAE3784B7, 0x0D8C34C1, 0xE7D54CF3);
	QapFp_set_u32(ifc->Rg->X->a0,
    0x061A10BB, 0x519EB62F, 0xEB8D8C7E, 0x8C61EDB6,
    0xA4648BBB, 0x4898BF0D, 0x91EE4224, 0xC803FB2B);
	
	QapFp_set_u32(ifc->Rg->Y->a1,
    0x0EBB2B0E, 0x7C8B1526, 0x8F6D4456, 0xF5F38D37,
    0xB09006FF, 0xD739C957, 0x8A2D1AEC, 0x6B3ACE9B);
	QapFp_set_u32(ifc->Rg->Y->a0,
    0x021897A0, 0x6BAF9343, 0x9A90E096, 0x698C8223,
    0x29BD0AE6, 0xBDBE09BD, 0x19F0E078, 0x91CD2B9A);
#else 
  QapFp_set_u32(ifc->Rg->X->a1,
    0, 0, 0, 0,
    0, 0, 0x11B4D9, 0xB9078C9C);
  QapFp_set_u32(ifc->Rg->X->a0,
    0, 0, 0, 0,
    0, 0, 0x1183E4, 0x2BE63029);

  QapFp_set_u32(ifc->Rg->Y->a1,
    0, 0, 0, 0,
    0, 0, 0x44923, 0xE174AB60);
  QapFp_set_u32(ifc->Rg->Y->a0,
    0, 0, 0, 0,
    0, 0, 0x10B2AC, 0x7252D760);
#endif // TINY

	// Create representations of zero, since we use it a lot
	QapFp_wp_init(ifc->Lzero);
	QapFp_copy(ifc->Lzero->X, ifc->Lg->X);
  QapFp_set_ui(ifc->Lzero->Y, 0); // Should be point at infinity in our internal representation

	QapFp2_wp_init(ifc->Rzero);
	QapFp2_copy(ifc->Rzero->X, ifc->Rg->X);
	QapFp2_set_ui(ifc->Rzero->Y, 0);
}

