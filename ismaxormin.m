
function out = ismaxormin(x)
    out = zeros(size(x));
    [~,maxi] = max(x);
    [~,mini] = min(x);
    out(maxi) = 1;
    out(mini) = 1;
end