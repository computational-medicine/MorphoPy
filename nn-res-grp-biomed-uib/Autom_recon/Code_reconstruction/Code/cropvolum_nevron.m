function [y,t1,t2,t3,t4,t5,t6]=cropvolum_nevron(stack)

[z1,z2,z3] = size(stack);
% t1
for i=1:z2
    t1=i;
    w1=find(stack(:,i,:)==1);
    if ~isempty(w1>=1)
        break
    end
end
%if t1>=2
%    t1=t1-1;
%end

% t2
for i=1:z2
    t2=z2-i+1;
    w2=find(stack(:,z2-i+1,:)==1);
    if ~isempty(w2>=1)
        break
    end
end
%if t2<=z2-1
%    t2=t2+1;
%end

% t3
for i=1:z1
    t3=i;
    w3=find(stack(i,:,:)==1);
    if ~isempty(w3>=1)
        break
    end
end
%if t3>=2
%    t3=t3-1;
%end

% t4
for i=1:z1
    t4=z1-i+1;
    w4=find(stack(z1-i+1,:,:)==1);
    if ~isempty(w4>=1)
        break
    end
end
%if t4<=z1-1
%    t4=t4+1;
%end

% t5
for j=1:z3
    t5=j;
    w5=find(stack(:,:,j)==1);
    if ~isempty(w5>=1)
        break
    end
end

% t6
for j=1:z3
    t6=z3-j+1;
    w6=find(stack(:,:,z3-j+1)==1);
    if ~isempty(w6>=1)
        break
    end
end

y=stack(t3:t4,t1:t2,t5:t6);