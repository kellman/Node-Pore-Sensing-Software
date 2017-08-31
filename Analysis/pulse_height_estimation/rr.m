function amp = rr(filt,y)


    cvx_begin quiet
        variable x
        minimize norm(filt*x-y,1)
        
    cvx_end
    amp = x;

end