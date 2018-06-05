function call = FDM_4(s0,k,vol,r,t,nprice,ntime,w)

call_1 = FDM_3(s0,k,vol,r,t,nprice,ntime,w);
call_2 = FDM_3(s0,k,vol,r,t,2*nprice,4*ntime,w);
call = (4*call_2 - call_1)/3;

end