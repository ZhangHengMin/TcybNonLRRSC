function B = Selection( A , k )


% B = A ;
% fprintf('K=%d,----------------------------\n',k) ;
if k == 0
    B = A ;
else
    [a, ind] = sort(A,'descend') ;
    B = zeros( size(A) ) ;    
    if k < 1        
        % ��ʱ��k��ʾ��������ϵ���ı���
%         fprintf('K=%d,----------------------------\n',k) ;

        k = round( size(A,1) * k ) ;
    end
    for i = 1 : size(A,2)
        B( ind(1:k,i) , i ) = A( ind(1:k,i) , i ) ;
    end
end



% ����0��ȫ����Ϊ1
% B(B>0) = 1 ;

% ����Ĳ�����YaleB����������Ч������Hopkins���½�Ч��

% for i = 1 : size(B,2)
%    B(:,i) = B(:,i) / max(B(:,i)) ;    
% end


