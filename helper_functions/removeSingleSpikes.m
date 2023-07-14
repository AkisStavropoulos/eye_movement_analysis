function x_clean = removeSingleSpikes(x,base_x,peak_x)

% applied on rectangular pulse function where you want to remove
% single-index rectangles

if all(size(x))>1; error('x cannot be a matrix'); end

x_clean = x;

for i = 1:numel(x)
    
    if i~=1
        
        if x(i)==peak_x && x(i-1)==base_x && x(i+1)==base_x
           
            x_clean(i) = base_x;
            
        elseif x(i)==base_x && x(i-1)==peak_x && x(i+1)==peak_x
            
            x_clean(i) = peak_x;
            
        end
        
    end    
    
end