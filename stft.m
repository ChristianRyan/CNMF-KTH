function stft = stft(frequence,M,R,N,w)  
number=ceil(N/R)+6;
input=zeros(number*R,1);
input(3*R+1:3*R+N)=frequence;
stft = [];
index=0;
%len=0;
while index<= length(input)-M
  fft_frame = fft(input(index+1:index+M).*w,M,1);
  stft = [stft [0;fft_frame(2:M/2+1,1)]];
  %stft = [stft fft_frame(1:M/2+1,1)];
  index=index+R;
  %len=len+M/2+1;
end
end