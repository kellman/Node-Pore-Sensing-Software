function amp = ls(filt,y)


    cvx_begin quiet
        variable x
        minimize norm(filt*x-y,2)
        
    cvx_end
    amp = x;

end