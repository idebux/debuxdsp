//---------------------------------------------------------------------------

#include <vcl.h>
#include <math.h>
#pragma hdrstop

#include "ufft.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)


union u_complex_16
{
  DWORD v;
  struct
  {
    short r;
    short i;
  };
};

typedef union u_complex_16 complex_16_t, *p_complex_16_t;

DWORD complex_16_construct ( short r, short i )
{
  complex_16_t ret = {0};
  ret.r = r;
  ret.i = i;
  return ret.v;
}
//---------------------------------------------------------------------------

DWORD complex_16_abs ( DWORD v )
{
  complex_16_t C = { v };
  return (DWORD) sqrt ( C.r*C.r + C.i*C.i );
}
//---------------------------------------------------------------------------

void complex_16_array_construct ( short r[],
                                   DWORD C[],
                                   size_t size )
{
  for ( size_t i = 0; i < size; i++ )
  {
    C[i] = complex_16_construct ( r[i], 0 );
  }
}
//---------------------------------------------------------------------------

void complex_16_express ( DWORD v,
                          short &r,
                          short &i )
{
  complex_16_t C = { v };
  r = C.r;
  i = C.i;
}
//---------------------------------------------------------------------------


DWORD complex_mul ( DWORD v1, DWORD v2 )
{
  complex_16_t V1 = { v1 },
               V2 = { v2 };
  complex_16_t R = { 0 };
  int r = ( int ) V1.r * V2.r - ( int ) V1.i * V2.i;
  int i = ( int ) V1.r * V2.i + ( int ) V1.i * V2.r;
  R.r = ( short )( r >> 15 );
  R.i = ( short )( i >> 15 );
  return R.v;
}
//---------------------------------------------------------------------------

DWORD complex_add ( DWORD v1, DWORD v2 )
{
  complex_16_t V1 = { v1 },
               V2 = { v2 };
  complex_16_t R = { 0 };
  R.r = V1.r + V2.r;
  R.i = V1.i + V2.i;
  return R.v;
}
//---------------------------------------------------------------------------

DWORD complex_sub ( DWORD v1, DWORD v2 )
{
  complex_16_t V1 = { v1 },
               V2 = { v2 };
  complex_16_t R = { 0 };
  R.r = V1.r - V2.r;
  R.i = V1.i - V2.i;
  return R.v;
}
//---------------------------------------------------------------------------

DWORD binner_product_p_16 ( DWORD s1, DWORD s2,
                            DWORD c1, DWORD c2 )
{
  return complex_add ( complex_mul ( s1, c1 ),
                       complex_mul ( s2, c2 ) );
}
//---------------------------------------------------------------------------

DWORD binner_product_n_16 ( DWORD s1, DWORD s2,
                            DWORD c1, DWORD c2 )
{
  return complex_sub ( complex_mul ( s1, c1 ),
                       complex_mul ( s2, c2 ) );
}   
//---------------------------------------------------------------------------

void fft8_round ( DWORD s8[],
                  DWORD c8[],
                  DWORD r8[] )
{
  r8[0] = binner_product_p_16 ( s8[0], s8[1], c8[0], c8[1] );
  r8[1] = binner_product_p_16 ( s8[2], s8[3], c8[2], c8[3] );
  r8[2] = binner_product_p_16 ( s8[4], s8[5], c8[4], c8[5] );
  r8[3] = binner_product_p_16 ( s8[6], s8[7], c8[6], c8[7] );
  r8[4] = binner_product_n_16 ( s8[0], s8[1], c8[0], c8[1] );
  r8[5] = binner_product_n_16 ( s8[2], s8[3], c8[2], c8[3] );
  r8[6] = binner_product_n_16 ( s8[4], s8[5], c8[4], c8[5] );
  r8[7] = binner_product_n_16 ( s8[6], s8[7], c8[6], c8[7] );
}
//---------------------------------------------------------------------------

void fft_round ( DWORD S[], DWORD C[],
                 DWORD R[], size_t N )
{
  size_t h = N >> 1;
  size_t l = 0,
         r = 0;
  for ( size_t i = 0; i < h; i++ )
  {
    l = i << 1;
    r = l  + 1;
    R [ i ] = binner_product_p_16 ( S[l], S[r], C[l], C[r] );
    R [i+h] = binner_product_n_16 ( S[l], S[r], C[l], C[r] );
  }
}                       
//---------------------------------------------------------------------------

int rounded ( double d )
{
   return ( int )( d * 2.0 + 1.0 ) >> 1;
}
//---------------------------------------------------------------------------

void gen_w8 ( DWORD w8[] )
{
  for ( int i = 0; i < 8; i++ )
  {
    w8[i] = complex_16_construct (
      ( short )( rounded ( cos ( 2 * M_PI / 8.0 * i ) * 32767.0 ) ),
      ( short )( rounded ( sin ( 2 * M_PI / 8.0 * i ) * 32767.0 ) ) );
  }
}
//---------------------------------------------------------------------------

void gen_w ( DWORD W[], size_t N )
{
  W[0] = 32767;
  for ( size_t i = 1; i < N; i++ )
  {
    W [i] = complex_16_construct (
      ( short )( rounded ( cos ( 2 * M_PI / N* i ) * 32767.999 ) ),
      ( short )( rounded ( sin ( 2 * M_PI / N * i ) * 32767.999 ) ) );
  }
}
//---------------------------------------------------------------------------

void gen_w512 ( DWORD w512[] )
{
  w512[0] = 32767;
  for ( int i = 1; i < 512; i++ )
  {
    w512[i] = complex_16_construct (
      ( short )( rounded ( cos ( 2 * M_PI * i / 512.0 ) * 32767.0 ) ),
      ( short )( rounded ( sin ( 2 * M_PI * i / 512.0 ) * 32767.0 ) ) );
  }
}
//---------------------------------------------------------------------------

void gen_fft8_coff_round1 ( DWORD w8[], DWORD f8cr1[] )
{
  f8cr1[0] = 32767;
  f8cr1[1] = w8[0];
  f8cr1[2] = 32767;
  f8cr1[3] = w8[0];
  f8cr1[4] = 32767;
  f8cr1[5] = w8[0];
  f8cr1[6] = 32767;
  f8cr1[7] = w8[0];
}
//---------------------------------------------------------------------------

void gen_fft8_coff_round2 ( DWORD w8[], DWORD f8cr2[] )
{
  f8cr2[0] = 32767;
  f8cr2[1] = w8[0];
  f8cr2[2] = 32767;
  f8cr2[3] = w8[0];
  f8cr2[4] = 32767;
  f8cr2[5] = w8[2];
  f8cr2[6] = 32767;
  f8cr2[7] = w8[2];
}
//---------------------------------------------------------------------------

void gen_fft8_coff_round3 ( DWORD w8[], DWORD f8cr3[] )
{
  f8cr3[0] = 32767;
  f8cr3[1] = w8[0];
  f8cr3[2] = 32767;
  f8cr3[3] = w8[1];
  f8cr3[4] = 32767;
  f8cr3[5] = w8[2];
  f8cr3[6] = 32767;
  f8cr3[7] = w8[3];
}
//---------------------------------------------------------------------------

void gen_fft512_coff_1 ( DWORD w512[], DWORD f512cr1[] )
{
  for ( int i = 0; i < 512; i++ )
  {
    f512cr1[i] = 32767;
  }
}
//---------------------------------------------------------------------------


void gen_fft512_coff_2 ( DWORD w512[], DWORD f512cr2[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr2 [ i*2 ] = 32767;
  }
  for ( int i = 0; i < 128; i++ )
  {
    f512cr2 [ i*2 + 1] = 32767;
    f512cr2 [ i*2 + 257] = w512[128];
  }
}
//---------------------------------------------------------------------------

void gen_fft512_coff_3 ( DWORD w512[], DWORD f512cr3[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr3 [ i*2 ] = 32767;
  }
  for ( int i = 0; i < 64; i++ )
  {
    f512cr3 [ i*2 + 1   ] = 32767;
    f512cr3 [ i*2 + 65 ] = w512[128];
    f512cr3 [ i*2 + 129 ] = 0x8000;
    f512cr3 [ i*2 + 193 ] = w512[384];
  }
}
//---------------------------------------------------------------------------

void gen_fft512_coff_4 ( DWORD w512[], DWORD f512cr4[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr4 [ i*2 ] = 32767;
  }
  /*
   *  分成8段，段号 j = 0-7， 每段值是 w512[ 64 * j ]
   */
  for ( int i = 0; i < 32; i++ )
  {
    f512cr4 [ i*2 +   1 ] = w512[0];
    f512cr4 [ i*2 +  65 ] = w512[64];
    f512cr4 [ i*2 + 129 ] = w512[128];
    f512cr4 [ i*2 + 193 ] = w512[192];
    f512cr4 [ i*2 + 257 ] = w512[256];
    f512cr4 [ i*2 + 321 ] = w512[320];
    f512cr4 [ i*2 + 385 ] = w512[384];
    f512cr4 [ i*2 + 449 ] = w512[448];
  }
}
//---------------------------------------------------------------------------

void gen_fft512_coff_5 ( DWORD w512[], DWORD f512cr5[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr5 [ i*2 ] = 32767;
  }
  for ( int i = 0; i < 16; i++ )
  {
    for ( int j = 0; j < 16; j++ )
    {
      f512cr5 [ i*2 + 1 + j * 32 ] = w512[j*32];
    }
  }
}
//---------------------------------------------------------------------------

void gen_fft512_coff_6 ( DWORD w512[], DWORD f512cr6[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr6 [ i*2 ] = 32767;
  }
  for ( int i = 0; i < 8; i++ )
  {
    for ( int j = 0; j < 32; j++ )
    {
      f512cr6 [ i*2 + 1 + j * 16 ] = w512[j*16];
    }
  }
}
//---------------------------------------------------------------------------

void gen_fft512_coff_7 ( DWORD w512[], DWORD f512cr7[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr7 [ i*2 ] = 32767;
  }
  for ( int i = 0; i < 4; i++ )
  {
    for ( int j = 0; j < 64; j++ )
    {
      f512cr7 [ i*2 + 1 + j * 8 ] = w512[j*8];
    }
  }
}
//---------------------------------------------------------------------------

void gen_fft512_coff_8 ( DWORD w512[], DWORD f512cr8[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr8 [ i*2 ] = 32767;
  }
  for ( int i = 0; i < 2; i++ )
  {
    for ( int j = 0; j < 128; j++ )
    {
      f512cr8 [ i*2 + 1 + j * 4 ] = w512[j*4];
    }
  }
}
//---------------------------------------------------------------------------

void gen_fft512_coff_9 ( DWORD w512[], DWORD f512cr9[] )
{
  for ( int i = 0; i < 256; i++ )
  {
    f512cr9 [ i*2 ] = 32767;
  }
  for ( int i = 0; i < 1; i++ )
  {
    for ( int j = 0; j < 256; j++ )
    {
      f512cr9 [ i*2 + 1 + j * 2 ] = w512[ j*2 ];
    }
  }
}
//---------------------------------------------------------------------------

size_t integer_log2 ( size_t n )
{
  size_t ret = 0;   
#pragma warn -8060
  while ( n = ( n >> 1 ) )  
#pragma warn +8060
  {
    ret ++;
  }
  return ret;
}
//---------------------------------------------------------------------------

void init_fft_coff ( LPDWORD fc[], size_t N )
{
  LPDWORD W = ( LPDWORD ) malloc ( sizeof( DWORD ) * N );
  W[0] = 32767;
  for ( size_t i = 1; i < N; i++ )
  {
    W[i] = complex_16_construct (
      ( short )( rounded ( cos ( 2 * M_PI * i / ( double ) N ) * 32767.999 ) ),
      ( short )( rounded ( sin ( 2 * M_PI * i / ( double ) N ) * 32767.999 ) ) );
  }
  size_t H = N >> 1;
  size_t J = 0;
  size_t n = N;
  LPDWORD fcr = NULL;
  int r = 0;
#pragma warn -8060
  while ( n = ( n >> 1 ) )
#pragma warn +8060
  {
    J = H / n;
    fcr = fc[r++];
    // 偶数角标永远是 32767
    for ( size_t h = 0; h < H; h++ )
    {
      fcr [ h*2 ] = 32767;
    }
    /*
     *  对于奇数角标，分成 J 段，每段有 n 个元素
     *  每段内，元素值相同
     *  设某段段号 j = 0 → J - 1，则元素的值都是 W [ j * n ]
     *  每段的每个元素
     */
    for ( size_t i = 0; i < n; i++ )
    {
      for ( size_t j = 0; j < J; j++ )
      {
        fcr [ i*2 + 1 + 2 * n * j ] = W [  n * j ];
      }
    }
  }
  free ( W );
}

void fft_shuf_t ( DWORD ti[], DWORD to[], size_t N )
{

  int weight[32] = { 0 };
  weight [0] = 0;
  for ( int i = 1; i < 32; i++ )
  {
    weight [i] = N / ( 1 << i );
  }
  for ( size_t i = 0; i < N; i++ )
  {
    size_t v = i;
    size_t b = 1;
    size_t r = 0;
    while ( v )
    {
      if ( v & 1 )
      {
        r += weight[b];
      }
      v >>= 1;
      b ++;
    }
    to[r] = ti[i];
  }
}
//---------------------------------------------------------------------------

void fft_shuf_f ( DWORD fi[], DWORD fo[], size_t N )
{
  fo [ 0 ] = fi [ 0 ];
  for ( size_t i = 1; i < N; i++ )
  {
    fo [ i ] = fi [ N - i ];
  }
}
//---------------------------------------------------------------------------

void fft (   DWORD  T[],  DWORD F[],
             DWORD II[], DWORD IO[],
           LPDWORD  C[], size_t N )
{
  LPDWORD ii = II;
  LPDWORD io = IO;
  LPDWORD it = NULL;
  size_t n = N;
  fft_shuf_t ( T, io, N );
  int r = 0;
#pragma warn -8060
  while ( n = ( n >> 1 ) )
#pragma warn +8060
  {
    it = io;
    io = ii;
    ii = it;
    fft_round ( ii, C[r++], io, N );
  }
  fft_shuf_f ( io, F, N );
}
//---------------------------------------------------------------------------

LPDWORD complex_16_array_create ( size_t N )

{
  return ( LPDWORD ) malloc ( sizeof( DWORD ) * N );
}
//---------------------------------------------------------------------------

void complex_16_array_delete ( LPDWORD array )

{
  free ( array );
}
//---------------------------------------------------------------------------

LPDWORD * fft_coff_create ( size_t N )
{
  size_t R = integer_log2 ( N );
  LPDWORD * ret = ( LPDWORD * ) malloc ( R * sizeof ( LPDWORD ) );
  for ( size_t r = 0; r < R; r++ )
  {
    ret[r] = complex_16_array_create ( N );
  }
  init_fft_coff ( ret, N );
  return ret;
}
//---------------------------------------------------------------------------

void fft_coff_delete ( LPDWORD * coff, size_t N )
{
  size_t R = integer_log2 ( N );
  for ( size_t r = 0; r < R; r++ )
  {
    complex_16_array_delete ( coff[r] );
  }
  free ( coff );
}
//---------------------------------------------------------------------------

void s_fcomplex::add ( fcomplex_p c1,
                       fcomplex_p c2,
                       fcomplex_p r  )
{
  r->m_r = c1->m_r + c2->m_r;
  r->m_i = c1->m_i + c2->m_i;
}
//---------------------------------------------------------------------------

void s_fcomplex::sub ( fcomplex_p c1,
                       fcomplex_p c2,
                       fcomplex_p r  )
{
  r->m_r = c1->m_r - c2->m_r;
  r->m_i = c1->m_i - c2->m_i;
}
//---------------------------------------------------------------------------

void s_fcomplex::mul ( fcomplex_p c1,
                       fcomplex_p c2,
                       fcomplex_p r  )
{
  r->m_r = c1->m_r * c2->m_r - c1->m_i * c2->m_i;
  r->m_i = c1->m_r * c2->m_i + c2->m_r * c1->m_i;
}
//---------------------------------------------------------------------------

void s_fcomplex::mul2 ( fcomplex_p c1,
                        fcomplex_p c2,
                        fcomplex_p r  )
{
  BYTE xmm[16] = {0};
  memcpy ( xmm, c1, 8 );
  memcpy ( xmm + 8, c2, 8 );
__asm
{
  // xmm0 = xx xx c2i c2r
  movdqu xmm1, dqword ptr [c2]


  // xmm0 = c2i c2r c1i c1r 
  movdqu xmm0, dqword ptr xmm
  // xmm1 = c2i c2i c1i c1r
  pshufd xmm1, xmm0, 11110100b
  // xmm0 = c1i c1r c2r c2r
  pshufd xmm0, xmm0, 01001010b
  // xmm1 = c1i*c2i c2i*c1r c1i*c2r c1r*c2r
  mulps xmm1, xmm0
  // xmm0 = c1r*c2r c2i*c1r c1i*c2r c1r*c2r
  pshufd xmm0, xmm1, 00100100b
  // xmm2 =      xx c1i*c2r      xx      xx
  pshufd xmm2, xmm1, 00010000b
  // xmm0 = c1r*c2r-c1i*c2i 00 00 00
  subps  xmm0, xmm1
  // xmm2 = xx c2i*c1r+c1i*c2r xx xx 
  addps  xmm2, xmm1
  // xmm0 = 00 00 00 c1r*c2r-c1i*c2i
  psrldq xmm0, 12
  // xmm2 = 00 xx c2i*c1r+c1i*c2r xx
  psrldq xmm2, 4
  // xmm2 = 00 00 c2i*c1r+c1i*c2r 00
  pshufd xmm1, xmm2, 11110111b
  por xmm0, xmm1
  movdqu dqword ptr xmm, xmm0
}
  memcpy ( r, xmm, 8 );
/*
  r->m_r = c1->m_r * c2->m_r - c1->m_i * c2->m_i;
  r->m_i = c1->m_r * c2->m_i + c2->m_r * c1->m_i;*/
}
//---------------------------------------------------------------------------

void s_fcomplex::mul_add( fcomplex_p s, fcomplex_p c,
                          fcomplex_p r )
{
  fcomplex_t c1,
             c2;
  mul2 ( s,   c,   &c1 );
  mul ( s+1, c+1, &c2 );
  add ( &c1, &c2, r );
}
//---------------------------------------------------------------------------

void s_fcomplex::mul_add_sub2( fcomplex_p s,  fcomplex_p c,
                               fcomplex_p psum, fcomplex_p pdef )
{
  float xmm[4];
__asm
{
  pushad
  mov eax, s
  mov ebx, c
  mov ecx, psum
  mov edx, pdef
  movdqu xmm0, dqword ptr [eax]
  movdqu dqword ptr xmm, xmm0   // debux
                                         movdqu xmm1, dqword ptr [ebx]   
  movdqu dqword ptr xmm, xmm1   // debux
                        movdqa xmm4,xmm0
  shufps xmm0,xmm1,01000100b
                                         shufps xmm4,xmm1,11101110b
  pshufd xmm1,xmm0, 11110100b            
                                         pshufd xmm5,xmm4, 11110100b
  pshufd xmm0,xmm0, 01001010b            
                                         pshufd xmm4,xmm4, 01001010b  
  mulps xmm0,xmm1                        
                                         mulps xmm4,xmm5   
  pxor xmm1,xmm1                         
                                         pxor xmm5,xmm5
  pxor xmm2,xmm2                                        
                                         pxor xmm6,xmm6
  shufps xmm1,xmm0,00100100b             
                                         shufps xmm5,xmm4,00100100b
  shufps xmm0,xmm2,11100111b             
                                         shufps xmm4,xmm6,11100111b
  pshufd xmm3,xmm0,11011110b             
                                         pshufd xmm7,xmm4,11011110b
  addps xmm1,xmm3                        
                                         addps xmm5,xmm7
  pshufd xmm0,xmm0,00101010b             
                                         pshufd xmm4,xmm4,00101010b
  subps xmm1,xmm0                        
                                         subps xmm5,xmm4
  pshufd xmm1,xmm1,01001011b
                                         pshufd xmm5,xmm5,01001011b
                        movdqa xmm0,xmm1
  subps xmm1,xmm5
                                         addps xmm0,xmm5
  movq qword ptr [edx], qword ptr xmm1
                                         movq qword ptr [ecx], qword ptr xmm0
  popad
}

}
//---------------------------------------------------------------------------

void s_fcomplex::mul_add_sub( fcomplex_p s, fcomplex_p c,
                              fcomplex_p psum, fcomplex_p pdef )
{
  fcomplex_t c1,
             c2;
  mul2 ( s,   c,   &c1 );
  mul2 ( s+1, c+1, &c2 );
  add ( &c1, &c2, psum );
  sub ( &c1, &c2, pdef );
}
//---------------------------------------------------------------------------

void s_fcomplex::mul_sub( fcomplex_p s, fcomplex_p c,
                          fcomplex_p r )
{
  fcomplex_t c1,
             c2;
  mul ( s,   c,   &c1 );
  mul ( s+1, c+1, &c2 );
  sub ( &c1, &c2, r );
}
//---------------------------------------------------------------------------

s_fcomplex::s_fcomplex ( float r, float i )
  : m_r ( r )
  , m_i ( i )
{
}
//---------------------------------------------------------------------------

s_fcomplex::s_fcomplex ( float r ) 
  : m_r ( r )
  , m_i ( 0 )
{
}
//---------------------------------------------------------------------------

s_fcomplex::s_fcomplex ( void )   
  : m_r ( 0 )
  , m_i ( 0 )
{
}
//---------------------------------------------------------------------------
                                                           
void s_fcomplex::set ( float r, float i )
{
  m_r = r;
  m_i = i;
}
//---------------------------------------------------------------------------

void s_fcomplex::set ( DWORD IQ )
{
  p_complex_16_t p = (p_complex_16_t)&IQ;
  m_r = (float)p->r;
  m_i = (float)p->i;
}
//---------------------------------------------------------------------------

void s_fcomplex::load_from_IQ16 ( fcomplex_p r, DWORD IQ[], size_t N )
{
  p_complex_16_t p = (p_complex_16_t)IQ;
  for ( size_t j = 0; j < N; j++ )
  {
    r[j].m_r = (float)p[j].r;
    r[j].m_i = (float)p[j].i;
  }
}
//---------------------------------------------------------------------------

void s_fcomplex::load_from_IQ16_ex ( fcomplex_p r, DWORD IQ[], size_t N )
{
__asm
{
  pushad
  mov ecx,N
  mov edi,r
  mov esi,IQ
  shr ecx,3
  jz end
  working:
  // {
  pxor xmm1,xmm1
                  pxor xmm2,xmm2
                                  pxor xmm5,xmm5
                                                  pxor xmm6,xmm6
  movdqu xmm0, dqword ptr [esi]
                                  movdqu xmm4, dqword ptr [esi+16]
  punpcklwd xmm1,xmm0
                  punpckhwd xmm2,xmm0
                                  punpcklwd xmm5,xmm4
                                                  punpckhwd xmm6,xmm4
  psrad xmm2,16
                  psrad xmm1,16
                                  psrad xmm5,16
                                                  psrad xmm6,16
  movdq2q mm0,xmm1
                  movdq2q mm2,xmm2
                                  movdq2q mm4,xmm5
                                                  movdq2q mm6,xmm6
  pshufd xmm1,xmm1, 01001110b
                  pshufd xmm2,xmm2, 01001110b
                                  pshufd xmm5,xmm5, 01001110b
                                                  pshufd xmm6,xmm6, 01001110b
  movdq2q mm1,xmm1
                  movdq2q mm3,xmm2
                                  movdq2q mm5,xmm5
                                                  movdq2q mm7,xmm6
  cvtpi2ps xmm0,mm0
  cvtpi2ps xmm4,mm1
                  cvtpi2ps xmm1,mm2
                  cvtpi2ps xmm5,mm3
                                  cvtpi2ps xmm2,mm4
                                  cvtpi2ps xmm6,mm5
                                                  cvtpi2ps xmm3,mm6
                                                  cvtpi2ps xmm7,mm7
  shufps xmm0,xmm4,01000100b
                  shufps xmm1,xmm5,01000100b
                                  shufps xmm2,xmm6,01000100b
                                                  shufps xmm3,xmm7,01000100b
  movdqu dqword ptr [edi], xmm0
                  movdqu dqword ptr [edi+16], xmm1
                                  movdqu dqword ptr [edi+32], xmm2
                                                  movdqu dqword ptr [edi+48], xmm3

  add esi,32        
  add edi,64        
  dec ecx
  jnz working       
  // }              
  end:              
  popad
  emms             
}
}
//---------------------------------------------------------------------------

fcomplex_p s_fcomplex::array_create ( size_t N )
{
  return ( fcomplex_p )malloc ( N * sizeof( fcomplex_t ) );
}
//---------------------------------------------------------------------------

void s_fcomplex::array_delete ( fcomplex_p array )
{
  free ( array );
}
//---------------------------------------------------------------------------

TFFT::TFFT ( size_t a_N )
  : m_N ( fix_N_if_error ( a_N ) )
  , m_c ( NULL )
  , m_ii ( NULL )
  , m_io ( NULL )
{
  set ( m_N );
}
//---------------------------------------------------------------------------

TFFT::~TFFT ()
{
  clr ();
}
//---------------------------------------------------------------------------

size_t TFFT::fix_N_if_error ( size_t a_N )
{
  return 1<<(integer_log2(a_N));
}
//---------------------------------------------------------------------------

void TFFT::set ( size_t N )
{
  clr ();
  m_ii = s_fcomplex::array_create ( m_N );
  m_io = s_fcomplex::array_create ( m_N );
   m_c = coff_create ( m_N );

  fcomplex_t w[8];
  memcpy ( w, m_c[0], 8 * sizeof( fcomplex_t ) );
  memcpy ( w, m_c[1], 8 * sizeof( fcomplex_t ) );
  memcpy ( w, m_c[2], 8 * sizeof( fcomplex_t ) );
}
//---------------------------------------------------------------------------

void TFFT::clr ( void )
{
  if ( m_ii )
  {
    s_fcomplex::array_delete ( m_ii );
  }
  if ( m_io )
  {
    s_fcomplex::array_delete ( m_io );
  }
  if ( m_c )
  {
    coff_delete ( m_c, m_N );
  }
}
//---------------------------------------------------------------------------

void TFFT::coff_init ( fcomplex_p c[], size_t N )
{
  fcomplex_p W = ( fcomplex_p ) malloc ( sizeof( s_fcomplex ) * N );
  for ( size_t i = 0; i < N; i++ )
  {
    W[i].m_r = (float) cos ( 2 * M_PI * i / ( double ) N );
    W[i].m_i = (float) sin ( 2 * M_PI * i / ( double ) N );
  }
  size_t H = N >> 1;
  size_t J = 0;
  size_t n = N;
  fcomplex_p fcr = NULL;
  int r = 0;
#pragma warn -8060
  while ( n = ( n >> 1 ) )
#pragma warn +8060
  {
    J = H / n;
    fcr = c[r++];
    // 偶数角标永远是 32767
    for ( size_t h = 0; h < H; h++ )
    {
      fcr [ h*2 ] = 1.0;
    }
    /*
     *  对于奇数角标，分成 J 段，每段有 n 个元素
     *  每段内，元素值相同
     *  设某段段号 j = 0 → J - 1，则元素的值都是 W [ j * n ]
     *  每段的每个元素
     */
    for ( size_t i = 0; i < n; i++ )
    {
      for ( size_t j = 0; j < J; j++ )
      {
        /*
         * 警告：结构体赋值！！！
         */
        fcr [ i*2 + 1 + 2 * n * j ] = W [  n * j ];
      }
    }
  }
  free ( W );
}
//---------------------------------------------------------------------------

fcomplex_p * TFFT::coff_create ( size_t N )
{
  size_t R = integer_log2 ( N );
  fcomplex_p * ret = ( fcomplex_p * ) malloc ( R * sizeof ( fcomplex_p ) );
  for ( size_t r = 0; r < R; r++ )
  {
    ret[r] = s_fcomplex::array_create ( N );
  }
  coff_init ( ret, N );
  return ret;
}
//---------------------------------------------------------------------------

void TFFT::coff_delete ( fcomplex_p *c, size_t N )
{
  size_t R = integer_log2 ( N );
  for ( size_t r = 0; r < R; r++ )
  {
    s_fcomplex::array_delete ( c[r] );
  }
  free ( c );
}
//---------------------------------------------------------------------------

void TFFT::round2 ( fcomplex_p S, fcomplex_p C, fcomplex_p R )
{
  size_t h = m_N >> 1;
  size_t l = 0;
  for ( size_t i = 0; i < h; i++ )
  {
    l = i << 1;
    s_fcomplex::mul_add_sub2 ( S+l, C+l, R+i, R+i+h );
  }
}
//---------------------------------------------------------------------------

void TFFT::round3 ( fcomplex_p S, fcomplex_p C, fcomplex_p R )
{
  size_t h = m_N >> 1;
  size_t l = 0;
  for ( size_t i = 0; i < h; i++ )
  {
    l = i << 1;
    s_fcomplex::mul_add_sub2 ( S+l, C+l, R+i, R+i+h );
  }
}
//---------------------------------------------------------------------------

void TFFT::round ( fcomplex_p S, fcomplex_p C, fcomplex_p R )
{
  size_t h = m_N >> 1;
  size_t l = 0;
  for ( size_t i = 0; i < h; i++ )
  {
    l = i << 1;
    s_fcomplex::mul_add ( S+l, C+l, R+i );
    s_fcomplex::mul_sub ( S+l, C+l, R+i+h );
  }
}
//---------------------------------------------------------------------------

void TFFT::shuf_t ( fcomplex_p ti, fcomplex_p to )
{
  int weight[32] = { 0 };
  weight [0] = 0;
  for ( int i = 1; i < 32; i++ )
  {
    weight [i] = m_N / ( 1 << i );
  }
  for ( size_t i = 0; i < m_N; i++ )
  {
    size_t v = i;
    size_t b = 1;
    size_t r = 0;
    while ( v )
    {
      if ( v & 1 )
      {
        r += weight[b];
      }
      v >>= 1;
      b ++;
    }
    to[r] = ti[i];
  }
}
//---------------------------------------------------------------------------

void TFFT::shuf_f ( fcomplex_p fi, fcomplex_p fo )
{
  fo [ 0 ] = fi [ 0 ];
  for ( size_t i = 1; i < m_N; i++ )
  {
    fo [ i ] = fi [ m_N - i ];
  }
}
//---------------------------------------------------------------------------

void TFFT::fft ( fcomplex_p t, fcomplex_p f )
{
  fcomplex_p ii = m_ii;
  fcomplex_p io = m_io;
  fcomplex_p it = NULL;
  size_t n = m_N;
  shuf_t ( t, io );
  int r = 0;
#pragma warn -8060
  while ( n = ( n >> 1 ) )
#pragma warn +8060
  {
    it = io;
    io = ii;
    ii = it;
    round ( ii, m_c[r++], io );
  }
  shuf_f ( io, f );
} 
//---------------------------------------------------------------------------

