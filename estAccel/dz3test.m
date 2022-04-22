function result = dz3test(vol)

s = [3.342604e-2 0.241125 0.450898 0.241125 3.342604e-2];
sd = [-9.186104e-2 -0.307610 0.00000 0.307610 9.186104e-2];
[x,y,z] = size(vol);
xConv = zeros(x,y,z);
for i = 1:y
    for j = 1:z
        xConv(:,i,j) = corrDn(vol(:,i,j),s','reflect1');
    end
end

xyConv = zeros(x,y,z);
for i = 1 : x
    for j = 1 : z
        xyConv(i,:,j) = corrDn(xConv(i,:,j),s,'reflect1');
    end
end

result = zeros(x,y,z);
for i = 1 : x
    for j = 1 : y
        result(i,j,:) = corrDn(xyConv(i,j,:),reshape(sd,[1 1 5]),'reflect1');
    end
end

return;

%%% Debug
vol=zeros([17,17,7]);
vol(:,9:17,:) = 1;
result3 = dz3(vol);
result2 = dx(vol(:,:,1));

end