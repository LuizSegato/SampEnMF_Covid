%This code calculates the SampEnMF for a given image

% [5] L. F. S. dos Santos, L. A. Neves, G. B. Rozendo, M. G. Ribeiro, M. Z. do Nascimento, T. A. A.
%Tosta, Multidimensional and fuzzy sample entropy (sampenmf) for quantifying h&e histological
%images of colorectal cancer, Computers in biology and medicine 103 (2018) 148–160.


%Test Healthy images
% ------------------------------------------

%General parameters
m_limit=4; %limit for parameter m (window size) 
k_limit=0.40; %limit for parameter k (tolerance constant)
k_start=0.06; %start of k values
k_increment=0.02; %increment for k values
limit_n=1602; %number of images this group

%Prepare the feature vector with (qtde_atr) attributes (variations of parameters m and k)
qtd_atr = floor(k_limit/k_increment)-floor(k_start/k_increment)+1;
SampEnMF=zeros(1,m_limit*qtd_atr);
matrix_SampEnMF=zeros(limit_n,m_limit*qtd_atr);
matrix_metrics=zeros(limit_n,m_limit*4+1);

%index image
for n=1:limit_n

  if isfile(strcat('../data/Healthy/Healthy (',num2str(n),').png')) 
      
    %Input image'
    image = imread(strcat('../data/Healthy/Healthy (',num2str(n),').png'));

    %attribute index
    i=1;
    
    %Sizes and amounts of sub-images
    sub_im_size=64; %sub-images size (64 x 64 pixels)
    N=size(image,1); %image height
    M=size(image,2); %image width
    sub_M=floor(M/sub_im_size); %amount of sub-images in the width of image
    sub_N=floor(N/sub_im_size);  %amount of sub-images in the height of image       

    %SampEnMF with multiscale approach
    %Sub-images preparation
    %Centering the initial position of the first sub-image considering the
    %amount of sub-images inside the image
    for m=1:m_limit
        for k=k_start:k_increment:k_limit
            %Necessary adjustment in height and width for scraps in a proportion of 25%
            if M-(sub_M*sub_im_size)>=floor(M/4) && N-(sub_N*sub_im_size)>=floor(N/4)
                im_x_from=floor((N-(sub_N*sub_im_size))/2);
                im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
                im_y_from=floor((M-(sub_M*sub_im_size))/2);
                im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
            %Necessary adjustment in width only for scraps in a proportion of 25%
            elseif M-(sub_M*sub_im_size)>=floor(M/4)
                im_x_from=1;
                im_x_to=sub_im_size;
                im_y_from=floor((M-(sub_M*sub_im_size))/2);
                im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
            %Necessary adjustment in height only for scraps in a proportion of 25%
            elseif N-(sub_N*sub_im_size)>=floor(N/4)
                im_x_from=floor((N-(sub_N*sub_im_size))/2);
                im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
                im_y_from=1;
                im_y_to=sub_im_size;
            %No necessary adjustment
            else
                im_x_from=1;
                im_x_to=sub_im_size;
                im_y_from=1;
                im_y_to=sub_im_size;
            end
        
            %Mean value of SampEnMF and other parameters
            SampEnMF_mean=0;
            qtd_black=0;
            entropy=0;
            sub_amount=(sub_M*sub_N); %amounts of sub-images inside image

            for im1=1:sub_M
                for im2=1:sub_N
                    %Calculate SampEnMF for each sub-image with parameters m and k
                    sub_image=image(im_x_from:im_x_to,im_y_from:im_y_to,:);
                    qtd_black=0;
                
                    for z=im_x_from:im_x_to
                        for l=im_y_from:im_y_to
                            if image(z,l,:) == 0
                                qtd_black=qtd_black+1;
                            end
                        end
                    end
                
                    %Consider backgroung above 90%
                    if qtd_black/(sub_im_size*sub_im_size)>=0.90 
                        sub_amount=sub_amount-1;
                    else
                        entropy=calcSampEnMF(sub_image,k,m);
                        SampEnMF_mean=SampEnMF_mean+double(entropy);
                    end
                    %Sub-image sliding in height
                    im_x_from=im_x_from+sub_im_size;
                    im_x_to=im_x_to+sub_im_size;
                end
                %Sub-image sliding in width
                im_y_from=im_y_from+sub_im_size;
                im_y_to=im_y_to+sub_im_size;
            
                %Re-centering the height case necessary 
                if N-(sub_N*sub_im_size)>=floor(N/4)
                    im_x_from=floor((N-(sub_N*sub_im_size))/2);
                    im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
                else
                    im_x_from=1;
                    im_x_to=sub_im_size;
                end
            end
            SampEnMF(i)=SampEnMF_mean/double(sub_amount);
            fprintf('Value of attribute %d (m = %d, e = %.2f) of SampEnMF = %.4f for healthy image %d\n',i,m,k,SampEnMF(i),n);
        
            %Next attribute
            i=i+1;
        end
    end

    %Calculate metrics on the SampEnMF curve (second composition of vector)
    %Example: m=4, 17 atrributes because the E (Maximum point scale) is the same for all m values
    metrics=zeros(1,m_limit*4+1);
    begin_curve=1;
    end_curve=5;
    begin_atrib=1;
    end_atrib=qtd_atr;
    k=k_start:k_increment:k_limit;

    for index_curves=1:m_limit
        metrics(1,begin_curve:end_curve)=curve_metrics(SampEnMF(1,begin_atrib:end_atrib),k_start,k_increment,k_limit);
    
        %Plot and save the SampEnMF curves for each scale m
        plot(k, SampEnMF(1, begin_atrib:end_atrib), 'b-s')
        xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
        ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
        title(strcat('m=',num2str(index_curves)))
    
        saveas(gcf, strcat('../results/curve_m',num2str(index_curves),'_healthy_(',num2str(n),').png'));
    
        aux=begin_curve;
        begin_curve=begin_curve+(end_curve-aux);
        end_curve=end_curve+(end_curve-aux);
        begin_atrib=begin_atrib+qtd_atr;
        end_atrib=end_atrib+qtd_atr;
    end    

    disp(metrics);
    
    %Matrix with all entropy and curves values
    matrix_SampEnMF(n,:)=SampEnMF;
    matrix_metrics(n,:)=metrics;
  end  
end

if m_limit>0
    xlswrite('../results/Healthy_attributes_matrix.csv',matrix_SampEnMF);
    xlswrite('../results/Healthy_metrics_matrix.csv',matrix_metrics);
end    

% -------------------------------------------


%Test Covid image
% ------------------------------------------

%General parameters
m_limit=4; %limit for parameter m (window size) 
k_limit=0.40; %limit for parameter k (tolerance constant)
k_start=0.06; %start of k values
k_increment=0.02; %increment for k values
limit_n=438; %number of images this group

%Prepare the feature vector with 24 attributes (variations of parameters m and k)
qtd_atr = floor(k_limit/k_increment)-floor(k_start/k_increment)+1;
SampEnMF=zeros(1,m_limit*qtd_atr); 
matrix_SampEnMF=zeros(limit_n,m_limit*qtd_atr);
matrix_metrics=zeros(limit_n,m_limit*4+1);

%index image
for n=1:limit_n

  if isfile(strcat('../data/Covid19/Covid (',num2str(n),').png'))   
      
    %Input image
    image = imread(strcat('../data/Covid19/Covid (',num2str(n),').png'));

    %attribute index
    i=1;
    
    %Sizes and amounts of sub-images
    sub_im_size=64; %sub-images size (64 x 64 pixels)
    N=size(image,1); %image height
    M=size(image,2); %image width
    sub_M=floor(M/sub_im_size); %amount of sub-images in the width of image
    sub_N=floor(N/sub_im_size);  %amount of sub-images in the height of image
       
    %SampEnMF with multiscale approach
    %Sub-images preparation
    %Centering the initial position of the first sub-image considering the
    %amount of sub-images inside the image
    for m=1:m_limit
        for k=k_start:k_increment:k_limit
            %Necessary adjustment in height and width for scraps in a proportion of 25%
            if M-(sub_M*sub_im_size)>=floor(M/4) && N-(sub_N*sub_im_size)>=floor(N/4)
                im_x_from=floor((N-(sub_N*sub_im_size))/2);
                im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
                im_y_from=floor((M-(sub_M*sub_im_size))/2);
                im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
            %Necessary adjustment in width only for scraps in a proportion of 25%
            elseif M-(sub_M*sub_im_size)>=floor(M/4)
                im_x_from=1;
                im_x_to=sub_im_size;
                im_y_from=floor((M-(sub_M*sub_im_size))/2);
                im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
            %Necessary adjustment in height only for scraps in a proportion of 25%
            elseif N-(sub_N*sub_im_size)>=floor(N/4)
                im_x_from=floor((N-(sub_N*sub_im_size))/2);
                im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
                im_y_from=1;
                im_y_to=sub_im_size;
            %No necessary adjustment
            else
                im_x_from=1;
                im_x_to=sub_im_size;
                im_y_from=1;
                im_y_to=sub_im_size;
            end
        
            %Mean value of SampEnMF and other parameters
            SampEnMF_mean=0;
            qtd_black=0;
            entropy=0;
            sub_amount=(sub_M*sub_N); %amounts of sub-images inside image
        
            for im1=1:sub_M
                for im2=1:sub_N
                    %Calculate SampEnMF for each sub-image with parameters m and k
                    sub_image=image(im_x_from:im_x_to,im_y_from:im_y_to,:);
                    qtd_black=0;
                
                    %If the sub-image is totally black (background) 
                    for z=im_x_from:im_x_to
                        for l=im_y_from:im_y_to
                            if image(z,l,:) == 0
                                qtd_black=qtd_black+1;
                            end
                        end
                    end
                
                    %Consider backgroung above 90%
                    if qtd_black/(sub_im_size*sub_im_size)>=0.90  
                        sub_amount=sub_amount-1;
                    else
                        entropy=calcSampEnMF(sub_image,k,m);
                        SampEnMF_mean=SampEnMF_mean+double(entropy);
                    end
                    %Sub-image sliding in height
                    im_x_from=im_x_from+sub_im_size;
                    im_x_to=im_x_to+sub_im_size;
                end
                %Sub-image sliding in width
                im_y_from=im_y_from+sub_im_size;
                im_y_to=im_y_to+sub_im_size;
            
                %Re-centering the height case necessary 
                if N-(sub_N*sub_im_size)>=floor(N/4)
                    im_x_from=floor((N-(sub_N*sub_im_size))/2);
                    im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
                else
                    im_x_from=1;
                    im_x_to=sub_im_size;
                end
            end
            SampEnMF(i)=SampEnMF_mean/double(sub_amount);
            fprintf('Value of attribute %d (m = %d, e = %.2f) of SampEnMF = %.4f for covid image %d\n',i,m,k,SampEnMF(i),n);
        
            %Next attribute
            i=i+1;
        end
    end

    %Calculate metrics on the SampEnMF curve (second composition of vector)
    %Example: m=4, 17 atrributes because the E (Maximum point scale) is the same for all m values
    metrics=zeros(1,m_limit*4+1);
    begin_curve=1;
    end_curve=5;
    begin_atrib=1;
    end_atrib=qtd_atr;
    k=k_start:k_increment:k_limit;

    for index_curves=1:m_limit
        metrics(1,begin_curve:end_curve)=curve_metrics(SampEnMF(1,begin_atrib:end_atrib),k_start,k_increment,k_limit);
    
        %Plot and save the SampEnMF curves for each scale m
        plot(k, SampEnMF(1, begin_atrib:end_atrib), 'b-s')
        xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
        ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
        title(strcat('m=',num2str(index_curves)))

        saveas(gcf, strcat('../results/curve_m',num2str(index_curves),'_covid_(',num2str(n),').png'));
    
        aux=begin_curve;
        begin_curve=begin_curve+(end_curve-aux);
        end_curve=end_curve+(end_curve-aux);
        begin_atrib=begin_atrib+qtd_atr;
        end_atrib=end_atrib+qtd_atr;
    end

    disp(metrics);
    
    %Matrix with all entropy and curve values
    matrix_SampEnMF(n,:)=SampEnMF;
    matrix_metrics(n,:)=metrics;
  end
end

if m_limit>0
    xlswrite('../results/Covid_attributes_matrix.csv',matrix_SampEnMF);
    xlswrite('../results/Covid_metrics_matrix.csv',matrix_metrics);
end
