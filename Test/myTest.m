%% 所有的单元测试都需要从matlab.unittest.TestCase继承
classdef myTest < matlab.unittest.TestCase
   
    %% 定义以Test为attribute的methods
    methods (Test)
        % 定义你自己的测试
        function testSingle(test) %function唯一的参数test是你的测试对象
            % Verifies single input case
                in        = single(10);             %输入
                expOut    = zeros(1,'single');      %期待的输出
                actualOut = foo(in);                %调用待测程序
                test.verifyEqual(actualOut,expOut); %比较实际输出与期待输出
        end
    end
end