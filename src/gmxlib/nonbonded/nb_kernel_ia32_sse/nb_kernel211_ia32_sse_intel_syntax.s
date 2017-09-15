;#
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as

.section .text



.globl nb_kernel211_ia32_sse
.globl _nb_kernel211_ia32_sse
nb_kernel211_ia32_sse:	
_nb_kernel211_ia32_sse:	
.equiv          nb211_p_nri,            8
.equiv          nb211_iinr,             12
.equiv          nb211_jindex,           16
.equiv          nb211_jjnr,             20
.equiv          nb211_shift,            24
.equiv          nb211_shiftvec,         28
.equiv          nb211_fshift,           32
.equiv          nb211_gid,              36
.equiv          nb211_pos,              40
.equiv          nb211_faction,          44
.equiv          nb211_charge,           48
.equiv          nb211_p_facel,          52
.equiv          nb211_argkrf,           56
.equiv          nb211_argcrf,           60
.equiv          nb211_Vc,               64
.equiv          nb211_type,             68
.equiv          nb211_p_ntype,          72
.equiv          nb211_vdwparam,         76
.equiv          nb211_Vvdw,             80
.equiv          nb211_p_tabscale,       84
.equiv          nb211_VFtab,            88
.equiv          nb211_invsqrta,         92
.equiv          nb211_dvda,             96
.equiv          nb211_p_gbtabscale,     100
.equiv          nb211_GBtab,            104
.equiv          nb211_p_nthreads,       108
.equiv          nb211_count,            112
.equiv          nb211_mtx,              116
.equiv          nb211_outeriter,        120
.equiv          nb211_inneriter,        124
.equiv          nb211_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb211_ixO,              0
.equiv          nb211_iyO,              16
.equiv          nb211_izO,              32
.equiv          nb211_ixH1,             48
.equiv          nb211_iyH1,             64
.equiv          nb211_izH1,             80
.equiv          nb211_ixH2,             96
.equiv          nb211_iyH2,             112
.equiv          nb211_izH2,             128
.equiv          nb211_iqO,              144
.equiv          nb211_iqH,              160
.equiv          nb211_dxO,              176
.equiv          nb211_dyO,              192
.equiv          nb211_dzO,              208
.equiv          nb211_dxH1,             224
.equiv          nb211_dyH1,             240
.equiv          nb211_dzH1,             256
.equiv          nb211_dxH2,             272
.equiv          nb211_dyH2,             288
.equiv          nb211_dzH2,             304
.equiv          nb211_qqO,              320
.equiv          nb211_qqH,              336
.equiv          nb211_c6,               352
.equiv          nb211_c12,              368
.equiv          nb211_six,              384
.equiv          nb211_twelve,           400
.equiv          nb211_vctot,            416
.equiv          nb211_Vvdwtot,          432
.equiv          nb211_fixO,             448
.equiv          nb211_fiyO,             464
.equiv          nb211_fizO,             480
.equiv          nb211_fixH1,            496
.equiv          nb211_fiyH1,            512
.equiv          nb211_fizH1,            528
.equiv          nb211_fixH2,            544
.equiv          nb211_fiyH2,            560
.equiv          nb211_fizH2,            576
.equiv          nb211_fjx,              592
.equiv          nb211_fjy,              608
.equiv          nb211_fjz,              624
.equiv          nb211_half,             640
.equiv          nb211_three,            656
.equiv          nb211_two,              672
.equiv          nb211_krf,              688
.equiv          nb211_crf,              704
.equiv          nb211_krsqO,            720
.equiv          nb211_krsqH1,           736
.equiv          nb211_krsqH2,           752
.equiv          nb211_is3,              768
.equiv          nb211_ii3,              772
.equiv          nb211_ntia,             776
.equiv          nb211_innerjjnr,        780
.equiv          nb211_innerk,           784
.equiv          nb211_n,                788
.equiv          nb211_nn1,              792
.equiv          nb211_nri,              796
.equiv          nb211_nouter,           800
.equiv          nb211_ninner,           804
.equiv          nb211_salign,           808
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 812		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb211_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb211_p_nri]
	mov ecx, [ecx]
	mov [esp + nb211_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb211_nouter], eax
	mov [esp + nb211_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb211_half], eax
	movss xmm1, [esp + nb211_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# 6.0
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# constant 12.0
	movaps [esp + nb211_half],  xmm1
	movaps [esp + nb211_two],  xmm2
	movaps [esp + nb211_three],  xmm3
	movaps [esp + nb211_six],  xmm4
	movaps [esp + nb211_twelve],  xmm5

	mov esi, [ebp + nb211_argkrf]
	mov edi, [ebp + nb211_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb211_krf], xmm5
	movaps [esp + nb211_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb211_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb211_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb211_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb211_iqO], xmm3
	movaps [esp + nb211_iqH], xmm4
	
	mov   edx, [ebp + nb211_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb211_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb211_ntia], ecx		

.nb211_threadloop:
        mov   esi, [ebp + nb211_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb211_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb211_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb211_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb211_n], eax
        mov [esp + nb211_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb211_outerstart
        jmp .nb211_end

.nb211_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb211_nouter]
	mov [esp + nb211_nouter], ebx

.nb211_outer:
	mov   eax, [ebp + nb211_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax +esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb211_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb211_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb211_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb211_pos]    ;# eax = base of pos[]  
	mov   [esp + nb211_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb211_ixO], xmm3
	movaps [esp + nb211_iyO], xmm4
	movaps [esp + nb211_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb211_ixH1], xmm0
	movaps [esp + nb211_iyH1], xmm1
	movaps [esp + nb211_izH1], xmm2
	movaps [esp + nb211_ixH2], xmm3
	movaps [esp + nb211_iyH2], xmm4
	movaps [esp + nb211_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb211_vctot], xmm4
	movaps [esp + nb211_Vvdwtot], xmm4
	movaps [esp + nb211_fixO], xmm4
	movaps [esp + nb211_fiyO], xmm4
	movaps [esp + nb211_fizO], xmm4
	movaps [esp + nb211_fixH1], xmm4
	movaps [esp + nb211_fiyH1], xmm4
	movaps [esp + nb211_fizH1], xmm4
	movaps [esp + nb211_fixH2], xmm4
	movaps [esp + nb211_fiyH2], xmm4
	movaps [esp + nb211_fizH2], xmm4
	
	mov   eax, [ebp + nb211_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb211_pos]
	mov   edi, [ebp + nb211_faction]	
	mov   eax, [ebp + nb211_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb211_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb211_ninner]
	mov   [esp + nb211_ninner], ecx
	add   edx, 0
	mov   [esp + nb211_innerk], edx    ;# number of innerloop atoms 
	jge   .nb211_unroll_loop
	jmp   .nb211_odd_inner
.nb211_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb211_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb211_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb211_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb211_iqO]
	mulps  xmm4, [esp + nb211_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb211_qqO], xmm3
	movaps  [esp + nb211_qqH], xmm4
	
	mov esi, [ebp + nb211_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb211_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb211_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb211_c6], xmm4
	movaps [esp + nb211_c12], xmm6

	mov esi, [ebp + nb211_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [esp + nb211_ixO]
	movaps xmm5, [esp + nb211_iyO]
	movaps xmm6, [esp + nb211_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb211_dxO], xmm4
	movaps [esp + nb211_dyO], xmm5
	movaps [esp + nb211_dzO], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb211_ixH1]
	movaps xmm5, [esp + nb211_iyH1]
	movaps xmm6, [esp + nb211_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb211_dxH1], xmm4
	movaps [esp + nb211_dyH1], xmm5
	movaps [esp + nb211_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb211_ixH2]
	movaps xmm4, [esp + nb211_iyH2]
	movaps xmm5, [esp + nb211_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [esp + nb211_dxH2], xmm3
	movaps [esp + nb211_dyH2], xmm4
	movaps [esp + nb211_dzH2], xmm5
	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + nb211_krf]	
	mulps  xmm1, [esp + nb211_krf]	
	mulps  xmm2, [esp + nb211_krf]	

	movaps [esp + nb211_krsqH2], xmm0
	movaps [esp + nb211_krsqH1], xmm1
	movaps [esp + nb211_krsqO], xmm2
	
	;# start with rsqO - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb211_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb211_half]
	movaps  xmm7, xmm4	;# rinvO in xmm7 
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb211_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb211_half]
	movaps  xmm6, xmm4	;# rinvH1 in xmm6 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb211_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb211_half]
	movaps  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb211_c6]
	mulps  xmm2, [esp + nb211_c12]
	movaps xmm3, xmm2
	subps  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addps  xmm3, [esp + nb211_Vvdwtot]
	mulps  xmm1, [esp + nb211_six]
	mulps  xmm2, [esp + nb211_twelve]
	subps  xmm2, xmm1	;# nb part of fs  

	movaps xmm0, xmm7
	movaps xmm1, [esp + nb211_krsqO]
	addps  xmm0, xmm1
	mulps  xmm1, [esp + nb211_two]
	subps  xmm0, [esp + nb211_crf] ;# xmm0=rinv+ krsq-crf 
	subps  xmm7, xmm1
	mulps  xmm0, [esp + nb211_qqO]
	mulps  xmm7, [esp + nb211_qqO]
	addps  xmm2, xmm7

	mulps  xmm4, xmm2	;# total fsO in xmm4 

	addps  xmm0, [esp + nb211_vctot]
	movaps [esp + nb211_Vvdwtot], xmm3
	movaps [esp + nb211_vctot], xmm0

	movaps xmm0, [esp + nb211_dxO]
	movaps xmm1, [esp + nb211_dyO]
	movaps xmm2, [esp + nb211_dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update O forces 
	movaps xmm3, [esp + nb211_fixO]
	movaps xmm4, [esp + nb211_fiyO]
	movaps xmm7, [esp + nb211_fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb211_fixO], xmm3
	movaps [esp + nb211_fiyO], xmm4
	movaps [esp + nb211_fizO], xmm7
	;# update j forces with water O 
	movaps [esp + nb211_fjx], xmm0
	movaps [esp + nb211_fjy], xmm1
	movaps [esp + nb211_fjz], xmm2

	;# H1 interactions 
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movaps  xmm7, xmm6
	movaps  xmm0, [esp + nb211_krsqH1]
	addps   xmm6, xmm0	;# xmm6=rinv+ krsq 
	mulps   xmm0, [esp + nb211_two]
	subps   xmm6, [esp + nb211_crf]
	subps   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulps   xmm6, [esp + nb211_qqH] ;# vcoul 
	mulps   xmm7, [esp + nb211_qqH]
	mulps  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addps  xmm6, [esp + nb211_vctot]

	movaps xmm0, [esp + nb211_dxH1]
	movaps xmm1, [esp + nb211_dyH1]
	movaps xmm2, [esp + nb211_dzH1]
	movaps [esp + nb211_vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H1 forces 
	movaps xmm3, [esp + nb211_fixH1]
	movaps xmm4, [esp + nb211_fiyH1]
	movaps xmm7, [esp + nb211_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb211_fixH1], xmm3
	movaps [esp + nb211_fiyH1], xmm4
	movaps [esp + nb211_fizH1], xmm7
	;# update j forces with water H1 
	addps  xmm0, [esp + nb211_fjx]
	addps  xmm1, [esp + nb211_fjy]
	addps  xmm2, [esp + nb211_fjz]
	movaps [esp + nb211_fjx], xmm0
	movaps [esp + nb211_fjy], xmm1
	movaps [esp + nb211_fjz], xmm2

	;# H2 interactions 
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movaps  xmm7, xmm5
	movaps  xmm0, [esp + nb211_krsqH2]
	addps   xmm5, xmm0	;# xmm5=rinv+ krsq 
	mulps   xmm0, [esp + nb211_two]
	subps   xmm5, [esp + nb211_crf]
	subps   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulps   xmm5, [esp + nb211_qqH] ;# vcoul 
	mulps   xmm7, [esp + nb211_qqH]
	mulps  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addps  xmm5, [esp + nb211_vctot]

	movaps xmm0, [esp + nb211_dxH2]
	movaps xmm1, [esp + nb211_dyH2]
	movaps xmm2, [esp + nb211_dzH2]
	movaps [esp + nb211_vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H2 forces 
	movaps xmm3, [esp + nb211_fixH2]
	movaps xmm4, [esp + nb211_fiyH2]
	movaps xmm7, [esp + nb211_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb211_fixH2], xmm3
	movaps [esp + nb211_fiyH2], xmm4
	movaps [esp + nb211_fizH2], xmm7

	mov edi, [ebp + nb211_faction]
	;# update j forces 
	addps xmm0, [esp + nb211_fjx]
	addps xmm1, [esp + nb211_fjy]
	addps xmm2, [esp + nb211_fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 136  ;# constant 10001000
	shufps xmm4, xmm7, 221  ;# constant 11011101			      
	;# xmm3 has fjx, xmm4 has fjy 
	subps xmm3, xmm0
	subps xmm4, xmm1
	;# unpack them back for storing 
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	;# finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 229  ;# constant 11100101
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 234  ;# constant 11101010
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 255  ;# constant 11111111
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb211_innerk],  4
	jl    .nb211_odd_inner
	jmp   .nb211_unroll_loop
.nb211_odd_inner:	
	add dword ptr [esp + nb211_innerk],  4
	jnz   .nb211_odd_loop
	jmp   .nb211_updateouterdata
.nb211_odd_loop:
	mov   edx, [esp + nb211_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb211_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb211_iqO]
	mov esi, [ebp + nb211_charge] 
	movhps xmm4, [esp + nb211_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb211_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov esi, [ebp + nb211_type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + nb211_vdwparam]
	shl ebx, 1	
	add ebx, [esp + nb211_ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [esp + nb211_c6], xmm6
	movaps [esp + nb211_c12], xmm7

	mov esi, [ebp + nb211_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb211_ixO]
	movss xmm4, [esp + nb211_iyO]
	movss xmm5, [esp + nb211_izO]
		
	movlps xmm6, [esp + nb211_ixH1]
	movlps xmm7, [esp + nb211_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb211_iyH1]
	movlps xmm7, [esp + nb211_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb211_izH1]
	movlps xmm7, [esp + nb211_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + nb211_dxO], xmm3
	movaps [esp + nb211_dyO], xmm4
	movaps [esp + nb211_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [esp + nb211_krf]
	movaps [esp + nb211_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb211_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb211_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000
		
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	mulss  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb211_c6]
	mulps  xmm2, [esp + nb211_c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb211_Vvdwtot]
	mulss  xmm1, [esp + nb211_six]
	mulss  xmm2, [esp + nb211_twelve]
	subss  xmm2, xmm1

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm3, [esp + nb211_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	mulps  xmm3, [esp + nb211_two]
	subps  xmm0, [esp + nb211_crf] ;# xmm0=rinv+ krsq-crf 
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [esp + nb211_qqO]	;# xmm0=vcoul 
	mulps  xmm1, [esp + nb211_qqO] 	;# xmm1=coul part of fs 

	addps xmm2, xmm1	;# total fs 
	
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addps  xmm0, [esp + nb211_vctot]
	movaps [esp + nb211_vctot], xmm0
	
	movaps xmm0, [esp + nb211_dxO]
	movaps xmm1, [esp + nb211_dyO]
	movaps xmm2, [esp + nb211_dzO]

	movaps [esp + nb211_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [esp + nb211_fixO]	
	movss  xmm4, [esp + nb211_fiyO]	
	movss  xmm5, [esp + nb211_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + nb211_fixO], xmm3	
	movss  [esp + nb211_fiyO], xmm4	
	movss  [esp + nb211_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [esp + nb211_fixH1]
	addss  xmm4, [esp + nb211_fiyH1]
	addss  xmm5, [esp + nb211_fizH1]
	movss  [esp + nb211_fixH1], xmm3	
	movss  [esp + nb211_fiyH1], xmm4	
	movss  [esp + nb211_fizH1], xmm5	;# updated the H1 force 

	mov edi, [ebp + nb211_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [esp + nb211_fixH2]
	addss  xmm4, [esp + nb211_fiyH2]
	addss  xmm5, [esp + nb211_fizH2]
	movss  [esp + nb211_fixH2], xmm3	
	movss  [esp + nb211_fiyH2], xmm4	
	movss  [esp + nb211_fizH2], xmm5	;# updated the H2 force 

	;# the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	;# x,y sum in xmm0 

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1 
	addss    xmm2, xmm1    ;# z sum in xmm2 
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec dword ptr [esp + nb211_innerk]
	jz    .nb211_updateouterdata
	jmp   .nb211_odd_loop
.nb211_updateouterdata:
	mov   ecx, [esp + nb211_ii3]
	mov   edi, [ebp + nb211_faction]
	mov   esi, [ebp + nb211_fshift]
	mov   edx, [esp + nb211_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb211_fixO]
	movaps xmm1, [esp + nb211_fiyO]
	movaps xmm2, [esp + nb211_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb211_fixH1]
	movaps xmm1, [esp + nb211_fiyH1]
	movaps xmm2, [esp + nb211_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb211_fixH2]
	movaps xmm1, [esp + nb211_fiyH2]
	movaps xmm2, [esp + nb211_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb211_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb211_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb211_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb211_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb211_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb211_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb211_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb211_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb211_n], esi
        jmp .nb211_outer
.nb211_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb211_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb211_end
        ;# non-zero, do one more workunit
        jmp   .nb211_threadloop
.nb211_end:
	emms

	mov eax, [esp + nb211_nouter]
	mov ebx, [esp + nb211_ninner]
	mov ecx, [ebp + nb211_outeriter]
	mov edx, [ebp + nb211_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb211_salign]
	add esp, eax
	add esp, 812
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



.globl nb_kernel211nf_ia32_sse
.globl _nb_kernel211nf_ia32_sse
nb_kernel211nf_ia32_sse:	
_nb_kernel211nf_ia32_sse:	
.equiv          nb211nf_p_nri,          8
.equiv          nb211nf_iinr,           12
.equiv          nb211nf_jindex,         16
.equiv          nb211nf_jjnr,           20
.equiv          nb211nf_shift,          24
.equiv          nb211nf_shiftvec,       28
.equiv          nb211nf_fshift,         32
.equiv          nb211nf_gid,            36
.equiv          nb211nf_pos,            40
.equiv          nb211nf_faction,        44
.equiv          nb211nf_charge,         48
.equiv          nb211nf_p_facel,        52
.equiv          nb211nf_argkrf,         56
.equiv          nb211nf_argcrf,         60
.equiv          nb211nf_Vc,             64
.equiv          nb211nf_type,           68
.equiv          nb211nf_p_ntype,        72
.equiv          nb211nf_vdwparam,       76
.equiv          nb211nf_Vvdw,           80
.equiv          nb211nf_p_tabscale,     84
.equiv          nb211nf_VFtab,          88
.equiv          nb211nf_invsqrta,       92
.equiv          nb211nf_dvda,           96
.equiv          nb211nf_p_gbtabscale,   100
.equiv          nb211nf_GBtab,          104
.equiv          nb211nf_p_nthreads,     108
.equiv          nb211nf_count,          112
.equiv          nb211nf_mtx,            116
.equiv          nb211nf_outeriter,      120
.equiv          nb211nf_inneriter,      124
.equiv          nb211nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb211nf_ixO,            0
.equiv          nb211nf_iyO,            16
.equiv          nb211nf_izO,            32
.equiv          nb211nf_ixH1,           48
.equiv          nb211nf_iyH1,           64
.equiv          nb211nf_izH1,           80
.equiv          nb211nf_ixH2,           96
.equiv          nb211nf_iyH2,           112
.equiv          nb211nf_izH2,           128
.equiv          nb211nf_iqO,            144
.equiv          nb211nf_iqH,            160
.equiv          nb211nf_qqO,            176
.equiv          nb211nf_qqH,            192
.equiv          nb211nf_c6,             208
.equiv          nb211nf_c12,            224
.equiv          nb211nf_vctot,          240
.equiv          nb211nf_Vvdwtot,        256
.equiv          nb211nf_half,           272
.equiv          nb211nf_three,          288
.equiv          nb211nf_krf,            304
.equiv          nb211nf_crf,            320
.equiv          nb211nf_krsqO,          336
.equiv          nb211nf_krsqH1,         352
.equiv          nb211nf_krsqH2,         368
.equiv          nb211nf_is3,            384
.equiv          nb211nf_ii3,            388
.equiv          nb211nf_ntia,           392
.equiv          nb211nf_innerjjnr,      396
.equiv          nb211nf_innerk,         400
.equiv          nb211nf_n,              404
.equiv          nb211nf_nn1,            408
.equiv          nb211nf_nri,            412
.equiv          nb211nf_nouter,         416
.equiv          nb211nf_ninner,         420
.equiv          nb211nf_salign,         424
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 428		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb211nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb211nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb211nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb211nf_nouter], eax
	mov [esp + nb211nf_ninner], eax


	mov esi, [ebp + nb211nf_argkrf]
	mov edi, [ebp + nb211nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb211nf_krf], xmm5
	movaps [esp + nb211nf_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb211nf_half], eax
	movss xmm1, [esp + nb211nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb211nf_half],  xmm1
	movaps [esp + nb211nf_three],  xmm3	

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb211nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb211nf_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb211nf_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb211nf_iqO], xmm3
	movaps [esp + nb211nf_iqH], xmm4
	
	mov   edx, [ebp + nb211nf_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov edi, [ebp + nb211nf_p_ntype]
	imul  ecx, [edi]      ;# ecx = ntia = 2*ntype*type[ii0] 
	mov   [esp + nb211nf_ntia], ecx		

.nb211nf_threadloop:
        mov   esi, [ebp + nb211nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb211nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb211nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb211nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb211nf_n], eax
        mov [esp + nb211nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb211nf_outerstart
        jmp .nb211nf_end

.nb211nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb211nf_nouter]
	mov [esp + nb211nf_nouter], ebx

.nb211nf_outer:
	mov   eax, [ebp + nb211nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb211nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb211nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb211nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb211nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb211nf_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb211nf_ixO], xmm3
	movaps [esp + nb211nf_iyO], xmm4
	movaps [esp + nb211nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb211nf_ixH1], xmm0
	movaps [esp + nb211nf_iyH1], xmm1
	movaps [esp + nb211nf_izH1], xmm2
	movaps [esp + nb211nf_ixH2], xmm3
	movaps [esp + nb211nf_iyH2], xmm4
	movaps [esp + nb211nf_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb211nf_vctot], xmm4
	movaps [esp + nb211nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb211nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb211nf_pos]
	mov   eax, [ebp + nb211nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb211nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb211nf_ninner]
	mov   [esp + nb211nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb211nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb211nf_unroll_loop
	jmp   .nb211nf_odd_inner
.nb211nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb211nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb211nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb211nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb211nf_iqO]
	mulps  xmm4, [esp + nb211nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb211nf_qqO], xmm3
	movaps  [esp + nb211nf_qqH], xmm4
	
	mov esi, [ebp + nb211nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb211nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb211nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb211nf_c6], xmm4
	movaps [esp + nb211nf_c12], xmm6

	mov esi, [ebp + nb211nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [esp + nb211nf_ixO]
	movaps xmm5, [esp + nb211nf_iyO]
	movaps xmm6, [esp + nb211nf_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2
	
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb211nf_ixH1]
	movaps xmm5, [esp + nb211nf_iyH1]
	movaps xmm6, [esp + nb211nf_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2
	
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb211nf_ixH2]
	movaps xmm4, [esp + nb211nf_iyH2]
	movaps xmm5, [esp + nb211nf_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + nb211nf_krf]	
	mulps  xmm1, [esp + nb211nf_krf]	
	mulps  xmm2, [esp + nb211nf_krf]	

	movaps [esp + nb211nf_krsqH2], xmm0
	movaps [esp + nb211nf_krsqH1], xmm1
	movaps [esp + nb211nf_krsqO], xmm2
	
	;# start with rsqO - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb211nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb211nf_half]
	movaps  xmm7, xmm4	;# rinvO in xmm7 
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb211nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb211nf_half]
	movaps  xmm6, xmm4	;# rinvH1 in xmm6 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb211nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb211nf_half]
	movaps  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb211nf_c6]
	mulps  xmm2, [esp + nb211nf_c12]
	movaps xmm3, xmm2
	subps  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addps  xmm3, [esp + nb211nf_Vvdwtot]

	movaps xmm0, xmm7
	movaps xmm1, [esp + nb211nf_krsqO]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb211nf_crf] ;# xmm0=rinv+ krsq-crf 
	subps  xmm7, xmm1
	mulps  xmm0, [esp + nb211nf_qqO]
	addps  xmm0, [esp + nb211nf_vctot]
	movaps [esp + nb211nf_Vvdwtot], xmm3
	movaps [esp + nb211nf_vctot], xmm0

	;# H1 interactions 
	movaps  xmm0, [esp + nb211nf_krsqH1]
	addps   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm6, [esp + nb211nf_crf]
	mulps   xmm6, [esp + nb211nf_qqH] ;# vcoul 
	addps  xmm6, [esp + nb211nf_vctot]
	movaps [esp + nb211nf_vctot], xmm6
	
	;# H2 interactions 
	movaps  xmm7, xmm5  ;# rinv 
	movaps  xmm0, [esp + nb211nf_krsqH2]
	addps   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subps   xmm5, [esp + nb211nf_crf]
	mulps   xmm5, [esp + nb211nf_qqH] ;# vcoul 
	addps   xmm5, [esp + nb211nf_vctot]
	movaps [esp + nb211nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb211nf_innerk],  4
	jl    .nb211nf_odd_inner
	jmp   .nb211nf_unroll_loop
.nb211nf_odd_inner:	
	add dword ptr [esp + nb211nf_innerk],  4
	jnz   .nb211nf_odd_loop
	jmp   .nb211nf_updateouterdata
.nb211nf_odd_loop:
	mov   edx, [esp + nb211nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb211nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb211nf_iqO]
	mov esi, [ebp + nb211nf_charge] 
	movhps xmm4, [esp + nb211nf_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb211nf_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov esi, [ebp + nb211nf_type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + nb211nf_vdwparam]
	shl ebx, 1	
	add ebx, [esp + nb211nf_ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [esp + nb211nf_c6], xmm6
	movaps [esp + nb211nf_c12], xmm7

	mov esi, [ebp + nb211nf_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb211nf_ixO]
	movss xmm4, [esp + nb211nf_iyO]
	movss xmm5, [esp + nb211nf_izO]
		
	movlps xmm6, [esp + nb211nf_ixH1]
	movlps xmm7, [esp + nb211nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb211nf_iyH1]
	movlps xmm7, [esp + nb211nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb211nf_izH1]
	movlps xmm7, [esp + nb211nf_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [esp + nb211nf_krf]
	movaps [esp + nb211nf_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb211nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb211nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000	

	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	mulss  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [esp + nb211nf_c6]
	mulps  xmm2, [esp + nb211nf_c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb211nf_Vvdwtot]
	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm3, [esp + nb211nf_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [esp + nb211nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [esp + nb211nf_qqO]	;# xmm0=vcoul 
	addps  xmm0, [esp + nb211nf_vctot]
	movaps [esp + nb211nf_vctot], xmm0
	movaps [esp + nb211nf_Vvdwtot], xmm5

	dec dword ptr [esp + nb211nf_innerk]
	jz    .nb211nf_updateouterdata
	jmp   .nb211nf_odd_loop
.nb211nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb211nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb211nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb211nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb211nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb211nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb211nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb211nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb211nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb211nf_n], esi
        jmp .nb211nf_outer
.nb211nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb211nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb211nf_end
        ;# non-zero, do one more workunit
        jmp   .nb211nf_threadloop
.nb211nf_end:
	emms

	mov eax, [esp + nb211nf_nouter]
	mov ebx, [esp + nb211nf_ninner]
	mov ecx, [ebp + nb211nf_outeriter]
	mov edx, [ebp + nb211nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb211nf_salign]
	add esp, eax
	add esp, 428
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

