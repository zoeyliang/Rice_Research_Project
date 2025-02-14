function result = L2err_toy(u1,u2,h)

result = sqrt(sum(h*(u1 - u2).^2));
end