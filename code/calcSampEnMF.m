%This function calculates the SampEnMF for a given color/grayscale image

% [5] L. F. S. dos Santos, L. A. Neves, G. B. Rozendo, M. G. Ribeiro, M. Z. do Nascimento, T. A. A.
%Tosta, Multidimensional and fuzzy sample entropy (sampenmf) for quantifying h&e histological
%images of colorectal cancer, Computers in biology and medicine 103 (2018) 148–160.

% Input:
% img - color image
% m - window size
% e - tolerance constant

% Output: 
% SampEnMF - a SampEnMF value (attribute)

function SampEnMF = calcSampEnMF(img,e,m)

%Convert uint8 to int16
img=im2int16(img);

N = size(img,1); %image height
M = size(img,2); %image width
C = size(img,3); %number of channels

tol=zeros(1,C);
r_tol=0;

if C==1
    r_tol = e*std2(img); %Calculate tolerance case grayscale image
else
    %Calculate tolerance case color image
    for color=1:C
        tol(color)=e*std2(img(:,:,color));
    end
end


Q = (N-m)*(M-m); %total number of windows inside image

Cmi=0;      %Average of the count over all windows compared to a fixed window (m x m pattern) (Equation 4 according to [5])
Cm1i=0;     %Average of the count over all windows compared to a fixed window (m+1 x m+1 pattern)
Cm=0;       %Similarity coefficient of all patterns analysed with size m × m (Equation 7 according to [5])
Cm1=0;      %similarity coefficient of all patterns analysed with size (m+1) × (m+1) (Equation 6 according to [5])
limit=180;   %Number of random windows (defined by [5])
rest=0;
%Linebase value to grayscale window
u0icinza=0;
u0jcinza=0; 

if mod(limit,9)>0
   rest=mod(limit,9); 
end

x=zeros(1,limit);
y=zeros(1,limit);

%Guaranteed random choice of windows from different regions of the image (9 regions)
im_x_from=1;
im_x_to=floor((N-m)/3);
im_y_from=1;
im_y_to=floor((M-m)/3);
vector_from=1;
vector_to=floor(limit/9);

%For number of windows >= 9 (at least one window for each region)
if limit>=9
    for z=1:3
        for w=1:3  
            rng('shuffle');
            x(1,vector_from:vector_to)=randi([im_x_from,im_x_to],1,floor(limit/9));
            y(1,vector_from:vector_to)=randi([im_y_from,im_y_to],1,floor(limit/9));
            vector_from=vector_from+floor(limit/9);
            vector_to=vector_to+floor(limit/9);
            im_x_from=im_x_from+floor((N-m)/3);
            im_x_to=im_x_to+floor((N-m)/3);
        end
        im_y_from=im_y_from+floor((M-m)/3);
        im_y_to=im_y_to+floor((M-m)/3);
        im_x_from=1;
        im_x_to=floor((N-m)/3);
    end    
end

%Adjusts for alocate the rest of windows in regions
im_x_from=1;
im_x_to=floor((N-m)/3);
im_y_from=1;
im_y_to=floor((M-m)/3);
position=vector_from;

for z=1:rest
    rng('shuffle');
    x(1,position)=randi([im_x_from,im_x_to]);
    y(1,position)=randi([im_y_from,im_y_to]);
    position=position+1;
    if mod((position-vector_from),3)==0
        im_x_from=1;
        im_x_to=floor((N-m)/3);
        im_y_from=im_y_from+floor((M-m)/3);
        im_y_to=im_y_to+floor((M-m)/3);
    else    
        im_x_from=im_x_from+floor((N-m)/3);
        im_x_to=im_x_to+floor((N-m)/3);
    end    
end

%Start quantification process for number of defined windows (limit)
%Parallel Computing Toolbox is used in this step (parfor)
parfor z=1:limit
    Cmi=0;Cm1i=0;
    %Linebase value to color window
    u0i=zeros(1,3); 
    u0j=zeros(1,3); 
    %Random window fixed for comparison with other windows in positions (i,j)
    j=x(z);
    i=y(z);
    
    %positions of other slided windows (i1,j1) -> rule: (i,j)~=(i1,j1)
    i1=i;
    j1=j+1;
         
    %Comparation process between windows
    while (i1<=N-m)
        while (j1<=M-m) 
            
            D=zeros(1,3);
            d=zeros(1,C);
            channel=zeros(1,3);
            
            %Case color image
            if C>1 
                
                  %Consider black background in a maximum of 2% of the window 
                  qtd_black_janela1 = 0;
                  qtd_black_janela2 = 0;
                
                  for k=1:m
                        for l=1:m 
                            if (img(i+k-1,j+l-1,:)) == -32768
                                qtd_black_janela1 = qtd_black_janela1 + 1;
                            end
                            if (img(i1+k-1,j1+l-1,:))  == -32768
                                qtd_black_janela2 = qtd_black_janela2 + 1;
                            end    
                        end       
                  end
                
                  if (qtd_black_janela1/(k*l)<=0.02) && (qtd_black_janela2/(k*l)<=0.02)
                      
                    for color=1:C
                        %Finding window baseline value for removal in Fuzzy approach (Equation 2 according [5])
                        u0i(color)=mean2(img(i:i+m-1,j:j+m-1,color)); 
                        u0j(color)=mean2(img(i1:i1+m-1,j1:j1+m-1,color));
                    end
                
                    %Similarity distance between windows
                    for k=1:m
                        for l=1:m
                            %SampEnMF with multidimensional approach
                            for color=1:C
                                %Equations 10 and 11 according [5]
                                d(color)=abs(abs(img(i+k-1,j+l-1,color)-u0i(color))-abs(img(i1+k-1,j1+l-1,color)-u0j(color)));
                            end  
                            d=double(d);
                            %Obtain maximum distance
                            [D(1,1),channel(1,1)]=max(d);
                            color_r=tol(channel(1,1));
                        end
                    end
                  
                    %SampEnMF with Fuzzy approach
                    D=double(D);
                    color_r=double(color_r);
                    if color_r==0
                        u=1;
                    else 
                        u=double(exp(-((double(D(1,1))^2)/(2*(color_r^2))))); %Equation 12 according [5]
                    end  
                        
                    Cmi=Cmi+u;
                    if (u~=0)
                        %Pattern m is similar, in this case the pattern m+1 is compared for the same window too
                        for l=1:m+1
                            %SampEnMF with multidimensional approach
                            for color=1:C
                                %Only the ultimate line of pattern m+1 is compared
                                d(color)=abs(abs(img(i+m,j+l-1,color)-u0i(color))-abs(img(i1+m,j1+l-1,color)-u0j(color)));
                            end
                            d=double(d);
                            %Obtain maximum distance
                            [D(1,2),channel(1,2)]=max(d);
                        end
                        for k=1:m+1
                            %SampEnMF with multidimensional approach
                            for color=1:C
                                %Only the ultimate column of pattern m+1 is compared
                                d(color)=abs(abs(img(i+k-1,j+m,color)-u0i(color))-abs(img(i1+k-1,j1+m,color)-u0j(color)));
                            end
                            d=double(d);
                            %Obtain maximum distance
                            [D(1,3),channel(1,3)]=max(d);
                        end
                    
                        [max_D,max_index]=max(D);
                        color_r=tol(channel(1,max_index));
                    
                        %SampEnMF with Fuzzy approach
                        max_D=double(max_D); 
                        color_r=double(color_r);
                        if color_r==0
                            u=1;
                        else 
                            u=double(exp(-((max_D^2)/(2*(color_r^2)))));
                        end
                      
                        Cm1i=Cm1i+u;
                    end
                  end  
                   
                %Case grayscale image
                else
                    
                  %Consider black background in a maximum of 2% of the window 
                  qtd_black_janela1 = 0;
                  qtd_black_janela2 = 0;
                  D1=0;
                
                  for k=1:m
                        for l=1:m 
                            if (img(i+k-1,j+l-1)) == -32768
                                qtd_black_janela1 = qtd_black_janela1 + 1;
                            end
                            if (img(i1+k-1,j1+l-1))  == -32768
                                qtd_black_janela2 = qtd_black_janela2 + 1;
                            end    
                        end       
                  end
                
                  if (qtd_black_janela1/(k*l)<=0.02) && (qtd_black_janela2/(k*l)<=0.02)
                    
                      %Finding window baseline value for removal in Fuzzy approach (Equation 2 according [5])
                      u0icinza=mean2(img(i:i+m-1,j:j+m-1));
                      u0jcinza=mean2(img(i1:i1+m-1,j1:j1+m-1));
                      for k=1:m
                          for l=1:m 
                              %SampEnMF with multidimensional approach
                              d=abs(abs(img(i+k-1,j+l-1)-u0icinza)-abs(img(i1+k-1,j1+l-1)-u0jcinza));   
                              d=double(d);
                              D1=max(d,D1);
                          end
                      end
                  
                      %SampEnMF with Fuzzy approach
                      D1=double(D1);
                      r=double(r_tol);
                      if r==0
                          u=1;
                      else    
                          u=double(exp(-((D1^2)/(2*(r^2)))));
                      end
                    
                      Cmi=Cmi+u;
                      if(u~=0)
                          %Pattern m is similar, in this case the pattern m+1 is compared for the same window too
                          for l=1:m+1
                              %SampEnMF with multidimensional approach
                              %Only the ultimate line of pattern m+1 is compared
                              d=abs(abs(img(i+m,j+l-1)-u0icinza)-abs(img(i1+m,j1+l-1)-u0jcinza));
                              d=double(d);
                              D1=max(d,D1);
                          end
                          for k=1:m+1
                              %SampEnMF with multidimensional approach
                              %Only the ultimate column of pattern m+1 is compared
                              d=abs(abs(img(i+k-1,j+m)-u0icinza)-abs(img(i1+k-1,j1+m)-u0jcinza));
                              d=double(d);
                              D1=max(d,D1);
                          end
                    
                          %SampEnMF with Fuzzy approach
                          D1=double(D1);
                          r=double(r_tol);
                          if r==0
                              u=1;
                          else 
                              u=double(exp(-((D1^2)/(2*(r^2)))));
                          end  
                      
                          Cm1i=Cm1i+u;
                      end
                  end  
            end
            %Slide window by columns
            j1=j1+1;
        end
        %Slide window by lines
        i1=i1+1;
        j1=1;
    end
    Cm=Cm+(double(Cmi)/(Q-1));
    Cm1=Cm1+(double(Cm1i)/(Q-1));
end
   
Cm=double(Cm)/Q;
Cm1=double(Cm1)/Q;

SampEnMF=abs(log(Cm1/Cm));