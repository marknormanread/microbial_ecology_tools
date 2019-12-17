"""


Mark N. Read, 2019
"""


import inspect
import plot_utils
import pandas as pd
import numpy as np
from os import makedirs
from shutil import rmtree
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import explained_variance_score


class PredictionResult:
    def __init__(self):
        self.grid: GridSearchCV = None
        self.predictors = None
        self.test_explained_variance = np.NAN
        self.test_permutation_pval = np.NAN  # P value under permutation testing, if performed.
        self.validation_explained_variance = np.NAN  # Training data, but assessed under cross validation
        self.validation_permutation_pval = np.NAN  # P value under permutation testing, if performed.
        self.y_train = []
        self.y_test = []
        self.y_validation_predictions = []
        self.y_test_predictions = []


def supervised_prediction_workflow(
        predictors: pd.DataFrame,
        response: list,
        available_data_as_test_proportion: float,  # Float [0, 1]
        pipeline: Pipeline,
        learning_mode: str,  # 'regression' or 'classification'
        parameters: dict,
        cv_mode=5,  # Number of cross validations to perform. Int, or 'loo' for leave one out.
        output_directory='',
        randomised_replicates=0,
        test_train_split_seed=3,
        save_hyperparamer_tuning_results=False,  # Should GridSearchCV acquire performance scores for all hyperparameter combos?
        figsize_mm=(80, 80),  # I work in mm, but matplotlib works in inches.
        # NOTE there appear to be problems with parallelism in sklearn. I've frequently found that the program
        # halts with no error messages if this is set to anything but "1". Googleing, it appears to be a somewhat
        # known problem relating to a multithreading backend in python.
        n_jobs=1,  # Number of cores to use. -1 uses all available.
        show_legend=True,
        validation_colour='#B05C0F',
        test_colour='#168C81',
        verbose_level=0  # Higher ints, more data to std out.
        ):
    """
    Example usage for pipeline and parameters:
    steps = [
        ('scaler', StandardScaler()),
        ('LASSO', Lasso())
    ]
    pipeline = Pipeline(steps)
    parameters = {
        'LASSO__alpha': np.logspace(-2, 2, 10)  # Note appending step name and '__' to actual parameter name.
    }
    """
    def perform_prediction(
            y,  # The response to be estimated. May be randomised.
            cv_mode,  # 'loo' or an integer specifying the number of folds.
            # Note that it makes little sense to be doing this for randomisation.
            save_hyperparamer_tuning_results
            ):
        result = PredictionResult()
        result.predictors = predictors

        # Stratify for classification, not for regression.
        stratify = None if learning_mode is 'regression' else y
        X_train, X_test, result.y_train, result.y_test = train_test_split(
            predictors, y, test_size=available_data_as_test_proportion,
            random_state=test_train_split_seed,
            stratify=stratify, shuffle=True
        )

        if cv_mode is 'loo':
            cv_mode = LeaveOneOut()

        # Find best hyperparameter values for the pipeline/estimator through cross validation.
        result.grid = GridSearchCV(  # This just sets up the experiment.
            pipeline, param_grid=parameters, cv=cv_mode,
            n_jobs=n_jobs,
            scoring='explained_variance',  # TODO this may not work for classification
            refit=True,  # Refit the estimator on the full training dataset using best parameters from CV.
            verbose=verbose_level,
            return_train_score=save_hyperparamer_tuning_results
        )
        result.grid.fit(X_train, result.y_train)  # This performs actual fitting

        # Evaluation against training data, but through cross validation.
        # IE, using estimator (model) parameters determined above, a model is fitted to the training portion of the
        # supplied data, and evaluated against the validation portion. All instances are used once in validation.
        result.y_validation_predictions = cross_val_predict(
            result.grid.best_estimator_,  # Can supply this directly, it retains best estimator hyperparameter values.
            X=X_train, y=result.y_train, cv=cv_mode, n_jobs=n_jobs
        )

        result.validation_explained_variance = \
            explained_variance_score(y_true=result.y_train, y_pred=result.y_validation_predictions)

        # Evaluate on independent test set (never used in model training)
        if len(result.y_test) > 0:
            result.y_test_predictions = result.grid.predict(X_test)

            result.test_explained_variance = \
                explained_variance_score(y_true=result.y_test, y_pred=result.y_test_predictions)
        return result

    def bootstrap_randomisation():
        bootstrap_results = pd.DataFrame(
            index=list(range(randomised_replicates)),
            columns=['validation_explained_variance', 'test_explained_variance']
        )
        for perm_i in list(range(randomised_replicates)):
            response_permutation = np.random.choice(response, size=len(response), replace=False)
            rand_result = perform_prediction(
                response_permutation, cv_mode=cv_mode,
                save_hyperparamer_tuning_results=False  # Nonsensical to do this under randomisation.
            )
            bootstrap_results.loc[perm_i, 'validation_explained_variance'] = rand_result.validation_explained_variance
            bootstrap_results.loc[perm_i, 'test_explained_variance'] = rand_result.test_explained_variance
        return bootstrap_results

    rmtree(output_directory, ignore_errors=True)  # Delete and make afresh.
    makedirs(output_directory, exist_ok=True)
    report = open(output_directory + '/report.txt', 'w')
    report_best_model = open(output_directory + '/best_model.txt', 'w')

    # Record the exact call made to this function.
    call_report = open(output_directory + '/call.txt', 'w')
    call_report.write('Prediction call:')
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    call_report.write('\n\tfunction name "%s"' % inspect.getframeinfo(frame)[2])
    for i in args:
        call_report.write("\n\t%s = %s" % (i, values[i]))
    call_report.close()

    true_result = perform_prediction(
        response, cv_mode=cv_mode,
        save_hyperparamer_tuning_results=save_hyperparamer_tuning_results
    )
    # This returns the score for each CV fold, and the mean across them.
    # It does not combine these data from across folds, as performed through cross_val_predict, below.
    if save_hyperparamer_tuning_results:
        pd.DataFrame(true_result.grid.cv_results_).to_csv(output_directory + '/cv_hyperparam_tuning_results.csv')

    report.write('\n\nBest model parameters:\n')
    report.write(str(true_result.grid.best_params_))
    report.write('\n\nBest estimator:\n')
    report.write(str(true_result.grid.best_estimator_))

    report_best_model.write('\n\nBest estimator:\n')
    report_best_model.write(str(true_result.grid.best_estimator_))

    report.write('\n\nValidation dataset performance (explained variance) = %3.2f'
                 % true_result.validation_explained_variance)
    report.write('\nTest dataset performance (explained variance) = %3.2f' % true_result.test_explained_variance)

    if randomised_replicates > 0:
        report.write('\n\nPerforming bootstrap randomisation.')
        bootstrap_results = bootstrap_randomisation()
        bootstrap_results.to_csv(output_directory + '/bootstrap_randomisation.csv',
                                 index_label='randomisation_replicate')
        p_val_resolution = 1. / randomised_replicates
        report.write('\nWith ' + str(randomised_replicates) + ' replicates, yields p val resolution of '
                     + str(p_val_resolution))

        rand_perf = bootstrap_results['validation_explained_variance']
        validation_beaten_by_random = rand_perf[rand_perf >= true_result.validation_explained_variance].count()
        rand_perf = bootstrap_results['test_explained_variance']
        test_beaten_by_random = rand_perf[rand_perf >= true_result.test_explained_variance].count()

        p_value_upper = (validation_beaten_by_random + 1) * p_val_resolution
        p_value_lower = validation_beaten_by_random * p_val_resolution
        true_result.validation_permutation_pval = p_value_upper
        report.write(
            '\nTrue VALIDATION data estimator performance = '
            + str(true_result.validation_explained_variance)
            + '\n  Randomisation beat this ' + str(validation_beaten_by_random) + ' times.'
            + '\n  P val in range ' + str(p_value_lower) + ' < p <= ' + str(p_value_upper)
        )
        p_value_upper = (test_beaten_by_random + 1) * p_val_resolution
        p_value_lower = test_beaten_by_random * p_val_resolution
        true_result.test_permutation_pval = p_value_upper
        report.write(
            '\n\nTrue TEST data estimator performance = '
            + str(true_result.test_explained_variance)
            + '\n  Randomisation beat this ' + str(test_beaten_by_random) + ' times.'
            + '\n  P val in range ' + str(p_value_lower) + ' < p <= ' + str(p_value_upper)
        )
        plot_utils.plot_cdfs(
            [bootstrap_results['validation_explained_variance'], bootstrap_results['test_explained_variance']],
            vertical_xs=[true_result.validation_explained_variance, true_result.test_explained_variance],
            names=['validation', 'test'],
            show_legend=show_legend,
            # colours=['royalblue', 'darkorange'],
            colours=[validation_colour, test_colour],
            xlabel='Explained variance',
            filename=output_directory + '/bootstrap_randomisation',
            fontsize=7,
            linewidth=1.0,
            figsize_mm=(50, 40)
        )

    # Graph performance, predicted vs actual response values.
    performance_df = pd.DataFrame({
        'Actual': list(true_result.y_train) + list(true_result.y_test),
        'Predicted': list(true_result.y_validation_predictions) + list(true_result.y_test_predictions),
        'Method': ['Validation'] * len(true_result.y_train) + ['Test'] * len(true_result.y_test)
    })
    f, ax = plt.subplots(
        figsize=[plot_utils.mm_to_in(figsize_mm[0]), plot_utils.mm_to_in(figsize_mm[1])],
        constrained_layout=True
    )
    ax = sns.scatterplot(
        data=performance_df,
        x='Actual', y='Predicted',
        hue='Method',
        palette={'Validation': validation_colour, 'Test': test_colour},
        linewidth=0,  # Disable marker outlines
        alpha=0.7,
        legend=show_legend,
        ax=ax  # Use axes with specified size.
    )
    # Draw line y=x; indicates exact relationship of actual to predicted response values.
    # Ideally everything lies on this line.
    plt.plot([np.min(response), np.max(response)], [np.min(response), np.max(response)], linewidth=2, c='k')
    plt.savefig(output_directory + '/actual_predicted.png', dpi=600)
    report.close()

    return true_result
