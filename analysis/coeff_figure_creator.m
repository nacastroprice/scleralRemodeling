% This will graph the coefficients over time

[r,c] = size(dtable);

coeff = zeros(r,11);
for i = 2:r
    for ii = 1:11
        coeff(i,ii) = dtable.coefficients{i,1}(ii);
    end
end


figure(1);
sgtitle("9122212: Ref inside inc")
subplot(5,2,1)
plot(dtable.dateTaken(:),coeff(:,1))
text(dtable.dateTaken(500),1,"a0")
ylim([0 4])
grid on
subplot(5,2,3)
plot(dtable.dateTaken(:),coeff(:,2))
text(dtable.dateTaken(500),0,"a1")
ylim([-12E-4 2E-4])
grid on
subplot(5,2,5)
plot(dtable.dateTaken(:),coeff(:,3))
text(dtable.dateTaken(500),0,"a2")
ylim([-12E-4 2E-4])
grid on
subplot(5,2,7)
plot(dtable.dateTaken(:),coeff(:,4))
text(dtable.dateTaken(500),0.6E-7,"a3")
grid on
ylim([-1.5E-7 1E-7])
subplot(5,2,9)
plot(dtable.dateTaken(:),coeff(:,5))
text(dtable.dateTaken(500),1E-7,"a4")
ylim([-4E-7 2E-7])
grid on
% figure(2);
subplot(5,2,2)
plot(dtable.dateTaken(:),coeff(:,6))
text(dtable.dateTaken(500),1,"b0")
ylim([0 4])
grid on
subplot(5,2,4)
plot(dtable.dateTaken(:),coeff(:,7))
text(dtable.dateTaken(500),1E-4,"b1")
ylim([-12E-4 2E-4])
grid on
subplot(5,2,6)
plot(dtable.dateTaken(:),coeff(:,8))
text(dtable.dateTaken(500),1.5E-4,"b2")
ylim([-12E-4 2E-4])
grid on
subplot(5,2,8)
plot(dtable.dateTaken(:),coeff(:,9))
text(dtable.dateTaken(500),0,"b3")
ylim([-1.5E-7 1E-7])
grid on
subplot(5,2,10)
plot(dtable.dateTaken(:),coeff(:,10))
text(dtable.dateTaken(500),0.4E-7,"b4")
ylim([-4E-7 2E-7])
grid on








