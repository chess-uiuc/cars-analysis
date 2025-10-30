classdef TestAddOne < matlab.unittest.TestCase
    methods (Test)
        function addsOneToScalar(testCase)
            act = add_one(41);
            testCase.verifyEqual(act, 42);
        end

        function vectorizedInput(testCase)
            act = add_one([0 1 2]);
            testCase.verifyEqual(act, [1 2 3]);
        end
    end
end
