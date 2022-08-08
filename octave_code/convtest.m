[s1,fs1] = audioread("G4.wav");
[s2,fs2] = audioread("gtrbody.wav");

sound(s1,fs1);
sound(s2,fs2);

t1 = 1/fs1:1/fs1:(1/fs1)*length(s1);        
t2 = 1/fs2:1/fs2:(1/fs2)*length(s2);
t1_ = 1/fs2:1/fs2:(1/fs1)*length(s1);
t2_ = 1/fs1:1/fs1:(1/fs2)*length(s2);

s1i = interp1(t1,s1,t1_);
s2i = interp1(t2,s2,t2_);


s3 = conv(s1i,s2);
s3norm = s3./max(s3);

sound(s3norm,fs2);
## sound(s1i,fs2);
## sound(s2i,fs1);
