function istft = istft(freq,M,R,N)
index=1;%0;
syn=[];
buf=zeros(M,1);
while index<=size(freq,2) %(size(freq,1)*size(freq,2))-M/2-1
    frame=freq(:,index);
    frame=[frame;flipud(conj(frame(2:end-1)))];
    ifft_frame = ifft(frame,M,1);
  if index==1
      syn = [syn ifft_frame];
  else
      ori=ifft_frame;
      %buf = [buf(R+1:end); zeros(R,1)];
      %buf = buf + ori;
      %syn=[syn;buf(1:R)];
      syn = [syn(1:end-3*R);syn(end-3*R+1:end)+ori(1:3*R);ori(end-R+1:end)];
   end
  index=index+1;
  %index=index+M/2+1;
end
istft = syn(3*R+1:3*R+N);
end
