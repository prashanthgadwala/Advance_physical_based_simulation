#include "fluid_solver.hpp"

#include <vislab/field/regular_field.hpp>
#include <vislab/graphics/colormap_texture.hpp>

namespace physsim
{
    FluidSolver::FluidSolver(const Eigen::AlignedBox2d& domain, const Eigen::Vector2i& resolution)
        : mDomain(domain)
        , mResolution(resolution)
        , mDensity(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mPressure(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mDivergence(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mVelocity_u(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mVelocity_v(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mAcceleration_u(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mAcceleration_v(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , densityTexture(std::make_shared<vislab::ColormapTexture>())
        , accuracy(1E-3)
        , iterations(1000)
        , stepSize(0.005 * std::sqrt((resolution.x() + resolution.y()) * 0.5))
        , pressureSolver(EPressureSolver::Iterative)
        , boundary(EBoundary::Closed)
    {
        mSpacing = (mDomain.max() - mDomain.min()).cwiseQuotient(mResolution.cast<double>());

        auto gridCentered = std::make_shared<vislab::RegularGrid2d>();
        gridCentered->setResolution(mResolution);
        gridCentered->setDomain(Eigen::AlignedBox2d(
            mDomain.min() + mSpacing * 0.5,
            mDomain.max() - mSpacing * 0.5));
        mDensity->setGrid(gridCentered);
        mDensity->setArray(std::make_shared<vislab::Array1f>());
        mDensity->getArray()->setSize(gridCentered->getResolution().prod());
        mPressure->setGrid(gridCentered);
        mPressure->setArray(std::make_shared<vislab::Array1f>());
        mPressure->getArray()->setSize(gridCentered->getResolution().prod());
        mDivergence->setGrid(gridCentered);
        mDivergence->setArray(std::make_shared<vislab::Array1f>());
        mDivergence->getArray()->setSize(gridCentered->getResolution().prod());

        auto gridStaggeredX = std::make_shared<vislab::RegularGrid2d>();
        gridStaggeredX->setResolution(mResolution + Eigen::Vector2i(1, 0));
        gridStaggeredX->setDomain(Eigen::AlignedBox2d(
            mDomain.min() + Eigen::Vector2d(0, mSpacing.y() * 0.5),
            mDomain.max() - Eigen::Vector2d(0, mSpacing.y() * 0.5)));
        mVelocity_u->setGrid(gridStaggeredX);
        mVelocity_u->setArray(std::make_shared<vislab::Array1f>());
        mVelocity_u->getArray()->setSize(gridStaggeredX->getResolution().prod());
        mAcceleration_u->setGrid(gridStaggeredX);
        mAcceleration_u->setArray(std::make_shared<vislab::Array1f>());
        mAcceleration_u->getArray()->setSize(gridStaggeredX->getResolution().prod());

        auto gridStaggeredY = std::make_shared<vislab::RegularGrid2d>();
        gridStaggeredY->setResolution(mResolution + Eigen::Vector2i(0, 1));
        gridStaggeredY->setDomain(Eigen::AlignedBox2d(
            mDomain.min() + Eigen::Vector2d(mSpacing.x() * 0.5, 0),
            mDomain.max() - Eigen::Vector2d(mSpacing.x() * 0.5, 0)));
        mVelocity_v->setGrid(gridStaggeredY);
        mVelocity_v->setArray(std::make_shared<vislab::Array1f>());
        mVelocity_v->getArray()->setSize(gridStaggeredY->getResolution().prod());
        mAcceleration_v->setGrid(gridStaggeredY);
        mAcceleration_v->setArray(std::make_shared<vislab::Array1f>());
        mAcceleration_v->getArray()->setSize(gridStaggeredY->getResolution().prod());

        mDensity_tmp    = std::shared_ptr<vislab::RegularSteadyScalarField2f>(mDensity->clone());
        mVelocity_u_tmp = std::shared_ptr<vislab::RegularSteadyScalarField2f>(mVelocity_u->clone());
        mVelocity_v_tmp = std::shared_ptr<vislab::RegularSteadyScalarField2f>(mVelocity_v->clone());

        densityTexture->scalarField               = mDensity;
        densityTexture->transferFunction.minValue = 0;
        densityTexture->transferFunction.maxValue = 1;
        densityTexture->transferFunction.values.clear();
        densityTexture->transferFunction.values.insert(std::make_pair(0., Eigen::Vector4d(84 / 255., 39 / 255., 143 / 255., 1)));
        densityTexture->transferFunction.values.insert(std::make_pair(0.25, Eigen::Vector4d(117 / 255., 107 / 255., 177 / 255., 1)));
        densityTexture->transferFunction.values.insert(std::make_pair(0.5, Eigen::Vector4d(158 / 255., 154 / 255., 200 / 255., 1)));
        densityTexture->transferFunction.values.insert(std::make_pair(0.75, Eigen::Vector4d(203 / 255., 201 / 255., 226 / 255., 1)));
        densityTexture->transferFunction.values.insert(std::make_pair(1., Eigen::Vector4d(242 / 255., 240 / 255., 247 / 255., 1)));
    }

    void FluidSolver::advance()
    {
        // apply source in density field
        applyDensitySource();

        // accumulate forces
        addBuoyancy();

        // advect the flow under application of the forces to the next time step. the flow becomes compressible.
        solveAdvection();

        // apply pressure-projection to make fluid incompressible
        solvePressure();

        // advect density in incompressible flow
        advectDensity();

        // for debugging, compute divergence
        computeDivergence();

        // reset forces
        mAcceleration_u->getArray()->setZero();
        mAcceleration_v->getArray()->setZero();
    }

    void FluidSolver::reset()
    {
        mDensity->getArray()->setZero();
        applyDensitySource();
        mPressure->getArray()->setZero();
        mDivergence->getArray()->setZero();
        mVelocity_u->getArray()->setZero();
        mVelocity_v->getArray()->setZero();
        mAcceleration_u->getArray()->setZero();
        mAcceleration_v->getArray()->setZero();

        // precompute sparse Cholesky factorization for Laplace operator on 2D grid
        Eigen::SparseMatrix<double> A(mResolution.prod(), mResolution.prod());
        buildLaplace2d(A);
        mPoissonFactorization.compute(A);
    }

    void FluidSolver::applyDensitySource()
    {
        Eigen::AlignedBox2d domain = mDensity->getDomain();
        Eigen::Vector2i resolution = mDensity->getGrid()->getResolution();
        Eigen::Vector2i low        = (sourceRegion.min() - domain.min()).cwiseQuotient(domain.max() - domain.min()).cwiseProduct(resolution.cast<double>() - Eigen::Vector2d::Ones()).cast<int>().cwiseMin(resolution - Eigen::Vector2i::Ones()).cwiseMax(Eigen::Vector2i::Zero());
        Eigen::Vector2i high       = (sourceRegion.max() - domain.min()).cwiseQuotient(domain.max() - domain.min()).cwiseProduct(resolution.cast<double>() - Eigen::Vector2d::Ones()).cast<int>().cwiseMin(resolution - Eigen::Vector2i::Ones()).cwiseMax(Eigen::Vector2i::Zero());

        for (int y = low.y(); y < high.y(); y++)
        {
            for (int x = low.x(); x < high.x(); x++)
            {
                float rho_ij = 1.f;
                mDensity->setVertexDataAt({ x, y }, rho_ij);
            }
        }
    }

    void FluidSolver::addBuoyancy()
    {
        Eigen::Vector2i resolution = mDensity->getGrid()->getResolution();
        float scaling              = 64.f / resolution.x();

        // approximate buoyancy simply from density
        for (int j = 1; j < resolution.y(); ++j)
        {
            for (int i = 0; i < resolution.x(); ++i)
            {
                float acc_v         = mAcceleration_v->getVertexDataAt({ i, j }).x();
                float density_below = mDensity->getVertexDataAt({ i, j - 1 }).x();
                float density_above = mDensity->getVertexDataAt({ i, j }).x();
                float density       = (density_below + density_above) / 2.f;
                acc_v += 0.1f * density * scaling;
                mAcceleration_v->setVertexDataAt({ i, j }, acc_v);
            }
        }
    }

    void FluidSolver::addAcceleration()
    {
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int j = 0; j < resu.y(); ++j)
        {
            for (int i = 0; i < resu.x(); ++i)
            {
                float vel_u = mVelocity_u->getVertexDataAt({ i, j }).x();
                float acc_u = mAcceleration_u->getVertexDataAt({ i, j }).x();
                vel_u += stepSize * acc_u;
                mVelocity_u->setVertexDataAt({ i, j }, vel_u);
            }
        }

        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int j = 0; j < resv.y(); ++j)
        {
            for (int i = 0; i < resv.x(); ++i)
            {
                float vel_v = mVelocity_v->getVertexDataAt({ i, j }).x();
                float acc_v = mAcceleration_v->getVertexDataAt({ i, j }).x();
                vel_v += stepSize * acc_v;
                mVelocity_v->setVertexDataAt({ i, j }, vel_v);
            }
        }
    }

    void FluidSolver::solveAdvection()
    {
        // accelerate the velocity field with the accumulated forces
        addAcceleration();

        // advect fluid -> the fluid becomes compressible
        advectVelocity();

        // set boundary conditions on velocity. it would be smarter if the advection procedure meets these automatically. In this form here, we lose energy.
        switch (boundary)
        {
        case EBoundary::Open:
            setNormalNeumann();
            break;
        case EBoundary::Closed:
            setNormalDirichlet();
            setTangentialNoSlip();
            break;
        }
    }

    void FluidSolver::solvePressure()
    {

        // compute divergence as input for the pressure solver
        computeDivergence();

        // compute pressure
        switch (pressureSolver)
        {
        case EPressureSolver::Iterative:
            solvePoissonIterative();
            break;
        case EPressureSolver::Linear:
            solvePoissonLinear();
            break;
        }

        // apply pressure to correct velocity
        correctVelocity();
    }

    void FluidSolver::setNormalNeumann()
    {
        // x-velocity
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int y = 0; y < resu.y(); ++y)
        {
            mVelocity_u->setVertexDataAt({ 0, y }, mVelocity_u->getVertexDataAt({ 1, y }));
            mVelocity_u->setVertexDataAt({ resu.x() - 1, y }, mVelocity_u->getVertexDataAt({ resu.x() - 2, y }));
        }

        // y-velocity
        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int x = 0; x < resv.x(); ++x)
        {
            mVelocity_v->setVertexDataAt({ x, 0 }, mVelocity_v->getVertexDataAt({ x, 1 }));
            mVelocity_v->setVertexDataAt({ x, resv.y() - 1 }, mVelocity_v->getVertexDataAt({ x, resv.y() - 2 }));
        }
    }

    void FluidSolver::setNormalDirichlet()
    {
        // x-velocity
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int y = 0; y < resu.y(); ++y)
        {
            mVelocity_u->setVertexDataAt({ 0, y }, 0);
            mVelocity_u->setVertexDataAt({ resu.x() - 1, y }, 0);
        }

        // y-velocity
        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int x = 0; x < resv.x(); ++x)
        {
            mVelocity_v->setVertexDataAt({ x, 0 }, 0);
            mVelocity_v->setVertexDataAt({ x, resv.y() - 1 }, 0);
        }
    }

    void FluidSolver::setTangentialNoSlip()
    {
        // x-velocity
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int x = 0; x < resu.x(); ++x)
        {
            mVelocity_u->setVertexDataAt({ x, 0 }, 0);
            mVelocity_u->setVertexDataAt({ x, resu.y() - 1 }, 0);
        }

        // y-velocity
        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int y = 0; y < resv.y(); ++y)
        {
            mVelocity_v->setVertexDataAt({ 0, y }, 0);
            mVelocity_v->setVertexDataAt({ resv.x() - 1, y }, 0);
        }
    }

    void FluidSolver::computeDivergence()
    {
        // calculate divergence
        Eigen::Vector2i res = mDivergence->getGrid()->getResolution();
        for (int y = 0; y < res.y(); ++y)
        {
            for (int x = 0; x < res.x(); ++x)
            {
                float xComponent = (mVelocity_u->getVertexDataAt({ x + 1, y }) - mVelocity_u->getVertexDataAt({ x, y })).x() / mSpacing.x();
                float yComponent = (mVelocity_v->getVertexDataAt({ x, y + 1 }) - mVelocity_v->getVertexDataAt({ x, y })).x() / mSpacing.y();
                float divergence = xComponent + yComponent;
                mDivergence->setVertexDataAt({ x, y }, divergence);
            }
        }
    }

    void FluidSolver::solvePoissonIterative()
    {
        float dx2      = mSpacing.prod();
        float residual = accuracy + 1; // initial residual
        float rho      = 1;

        for (int it = 0; residual > accuracy && it < iterations; ++it)
        {
            switch (boundary)
            {
            case EBoundary::Open:
                for (int j = 0; j < mResolution.y(); ++j)
                {
                    for (int i = 0; i < mResolution.x(); ++i)
                    {
                        float b = mDivergence->getVertexDataAt({ i, j }).x() / stepSize * rho; // right-hand
                        // TODO: update the pressure values
                        float p_left = (i > 0) ? mPressure->getVertexDataAt({ i - 1, j }).x() : 0;
                        float p_right = (i < mResolution.x() - 1) ? mPressure->getVertexDataAt({ i + 1, j }).x() : 0;
                        float p_bottom = (j > 0) ? mPressure->getVertexDataAt({ i, j - 1 }).x() : 0;
                        float p_top = (j < mResolution.y() - 1) ? mPressure->getVertexDataAt({ i, j + 1 }).x() : 0;
                        
                        float new_pressure = (p_left + p_right + p_bottom + p_top - dx2 * b) / 4.0f;
                        mPressure->setVertexDataAt({ i, j }, new_pressure);
                    }
                }
                break;
            case EBoundary::Closed:
                for (int j = 0; j < mResolution.y(); ++j)
                {
                    for (int i = 0; i < mResolution.x(); ++i)
                    {
                        float b = mDivergence->getVertexDataAt({ i, j }).x() / stepSize * rho; // right-hand
                        // TODO: update the pressure values
                        float p_left = (i > 0) ? mPressure->getVertexDataAt({ i - 1, j }).x() : mPressure->getVertexDataAt({ i, j }).x();
                        float p_right = (i < mResolution.x() - 1) ? mPressure->getVertexDataAt({ i + 1, j }).x() : mPressure->getVertexDataAt({ i, j }).x();
                        float p_bottom = (j > 0) ? mPressure->getVertexDataAt({ i, j - 1 }).x() : mPressure->getVertexDataAt({ i, j }).x();
                        float p_top = (j < mResolution.y() - 1) ? mPressure->getVertexDataAt({ i, j + 1 }).x() : mPressure->getVertexDataAt({ i, j }).x();
                        
                        int k = 0;
                        if (i > 0) k++;
                        if (i < mResolution.x() - 1) k++;
                        if (j > 0) k++;
                        if (j < mResolution.y() - 1) k++;
                        
                        float new_pressure = (p_left + p_right + p_bottom + p_top - dx2 * b) / k;
                        mPressure->setVertexDataAt({ i, j }, new_pressure);
                    }
                }
                break;
            }

            // Compute the new residual, i.e. the sum of the squares of the individual residuals (squared L2-norm)
            residual = 0;
            for (int j = 1; j < mResolution.y() - 1; ++j)
            {
                for (int i = 1; i < mResolution.x() - 1; ++i)
                {
                    float b = mDivergence->getVertexDataAt({ i, j }).x() / stepSize * rho; // right-hand
                    // TODO: compute the cell residual
                    float p_center = mPressure->getVertexDataAt({ i, j }).x();
                    
                    if (boundary == EBoundary::Open) {
                        // For open boundaries, use simple Laplacian with zero boundaries
                        float p_left = (i > 0) ? mPressure->getVertexDataAt({ i - 1, j }).x() : 0;
                        float p_right = (i < mResolution.x() - 1) ? mPressure->getVertexDataAt({ i + 1, j }).x() : 0;
                        float p_bottom = (j > 0) ? mPressure->getVertexDataAt({ i, j - 1 }).x() : 0;
                        float p_top = (j < mResolution.y() - 1) ? mPressure->getVertexDataAt({ i, j + 1 }).x() : 0;
                        float cell_residual = b - (-4*p_center + p_left + p_right + p_bottom + p_top) / dx2;
                        residual += cell_residual * cell_residual; // accumulate squared residual

                    } else {
                        // For closed boundaries, use Neumann-consistent Laplacian
                        float p_left = (i > 0) ? mPressure->getVertexDataAt({ i - 1, j }).x() : p_center;
                        float p_right = (i < mResolution.x() - 1) ? mPressure->getVertexDataAt({ i + 1, j }).x() : p_center;
                        float p_bottom = (j > 0) ? mPressure->getVertexDataAt({ i, j - 1 }).x() : p_center;
                        float p_top = (j < mResolution.y() - 1) ? mPressure->getVertexDataAt({ i, j + 1 }).x() : p_center;
                        
                        int k = 0;
                        if (i > 0) k++;
                        if (i < mResolution.x() - 1) k++;
                        if (j > 0) k++;
                        if (j < mResolution.y() - 1) k++;

                        float cell_residual = b - (-k * p_center + p_left + p_right + p_bottom + p_top) / dx2;
                        residual += cell_residual * cell_residual; // accumulate squared residual
                    }

                }
            }

            // Get the L2-norm of the residual
            residual = sqrt(residual);

            // We assume the accuracy is meant for the average L2-norm per grid cell
            residual /= (mResolution.x() - 2) * (mResolution.y() - 2);
        }
    }

    void FluidSolver::solvePoissonLinear()
    {
        float dx2 = mSpacing.prod();
        float rho = 1;

        // right hand side
        Eigen::VectorXd b(mResolution.x() * mResolution.y());
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                b(j * (std::size_t)mResolution.x() + i) = mDivergence->getVertexDataAt({ i, j }).x() / stepSize * rho * dx2;

        // solve sparse linear SPD system
        Eigen::VectorXd pout = mPoissonFactorization.solve(b);

        // copy result to output
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                mPressure->setVertexDataAt({ i, j }, pout(j * (std::size_t)mResolution.x() + i));
    }

    void FluidSolver::correctVelocity()
    {
        //     the staggered grid indices look like this
        //     ------------------------------
        //     |                            |
        //     |u_i-0.5,j      p_i,j        |u_i+0.5,j
        //     |                            |
        //     ------------------------------

        //     velocity u is actually stored at
        //     ------------------------------
        //     |                            |
        //     |u_i,j          p_i,j        |u_i+1,j
        //     |                            |
        //     ------------------------------
        //     hence, updating u_i,j needs p_i,j and p_i-1,j

        // Note: velocity u_{i+1/2} is practically stored at i+1, hence xV_{i} -= dt * (p_{i} - p_{i-1}) / dx
        for (int j = 1; j < mResolution.y() - 1; ++j)
            for (int i = 1; i < mResolution.x(); ++i)
            {
                // TODO: update u
                float p_current = mPressure->getVertexDataAt({ i, j }).x();
                float p_prev = mPressure->getVertexDataAt({ i - 1, j }).x();
                float u_current = mVelocity_u->getVertexDataAt({ i, j }).x();
                float rho = 1.0f;
                
                u_current -= stepSize * (p_current - p_prev) / mSpacing.x() / rho;
                mVelocity_u->setVertexDataAt({ i, j }, u_current);
            }

        // Same for velocity v_{i+1/2}.
        for (int j = 1; j < mResolution.y(); ++j)
            for (int i = 1; i < mResolution.x() - 1; ++i)
            {
                // TODO: update v
                float p_current = mPressure->getVertexDataAt({ i, j }).x();
                float p_prev = mPressure->getVertexDataAt({ i, j - 1 }).x();
                float v_current = mVelocity_v->getVertexDataAt({ i, j }).x();
                float rho = 1.0f;
                
                v_current -= stepSize * (p_current - p_prev) / mSpacing.y() / rho;
                mVelocity_v->setVertexDataAt({ i, j }, v_current);
            }
    }

    void FluidSolver::advectVelocity()
    {
        // Velocities live on the MAC grid
        // velocity_u resolution is: [0,m_res_x] x [0,m_res_y-1]
        // velocity_v resolution is: [0,m_res_x-1] x [0,m_res_y]

        //     the staggered grid indices look like this (with half offsets added already)
        //        v_i-1,j+1      v_i,j+1
        //     ------O-------------O--------
        //     |            |
        //     |            O uij  X pij
        //     |            |
        //     ------O-------------O--------
        //        v_i-1,j        v_i,j

        // Velocities (u), MAC grid
        for (int j = 0; j < mResolution.y(); ++j)
        {
            for (int i = 1; i < mResolution.x(); ++i)
            { // skip first and last row: those are determined by the boundary condition
                // TODO: Compute the velocity
                float u_current = mVelocity_u->getVertexDataAt({ i, j }).x();
                
                float v_left = (mVelocity_v->getVertexDataAt({ i - 1, j + 1 }).x() + 
                               mVelocity_v->getVertexDataAt({ i - 1, j }).x()) / 2.0f;
                float v_right = (mVelocity_v->getVertexDataAt({ i , j + 1 }).x() + 
                                mVelocity_v->getVertexDataAt({ i, j }).x()) / 2.0f;

                float last_x_velocity = u_current;
                float last_y_velocity = (v_left + v_right) / 2.0f;

                // TODO: Find the last position of the particle (in grid coordinates) using an Euler step
                float last_x = i - stepSize * last_x_velocity / mSpacing.x();
                float last_y = j - stepSize * last_y_velocity / mSpacing.y();

                // Make sure the coordinates are inside the boundaries
                if (last_x < 0)
                    last_x = 0;
                if (last_y < 0)
                    last_y = 0;
                if (last_x > mResolution.x() - 0)
                    last_x = mResolution.x() - 0;
                if (last_y > mResolution.y() - 1)
                    last_y = mResolution.y() - 1;

                // Determine corners for bilinear interpolation
                int x_low  = (int)last_x;
                int y_low  = (int)last_y;
                int x_high = std::min(x_low + 1, mResolution.x());
                int y_high = std::min(y_low + 1, mResolution.y() - 1);

                // Compute the interpolation weights
                float x_weight = last_x - x_low;
                float y_weight = last_y - y_low;

                // TODO: Bilinear interpolation
                float u00 = mVelocity_u->getVertexDataAt({ x_low, y_low }).x();
                float u10 = mVelocity_u->getVertexDataAt({ x_high, y_low }).x();
                float u01 = mVelocity_u->getVertexDataAt({ x_low, y_high }).x();
                float u11 = mVelocity_u->getVertexDataAt({ x_high, y_high }).x();
                
                float interpolated_u = (1 - x_weight) * (1 - y_weight) * u00 +
                                      x_weight * (1 - y_weight) * u10 +
                                      (1 - x_weight) * y_weight * u01 +
                                      x_weight * y_weight * u11;
                
                mVelocity_u_tmp->setVertexDataAt({ i, j }, interpolated_u);
            }
        }

        //     the staggered grid indices look like this (with half offsets added already)
        //     |              |
        //     O uij  X pij   | u_i+1,j
        //     |              |
        //     -------O--------
        //     |    v_i,j     |
        //     O uij-1        | u_i+1,j-1
        //     |              |

        // Velocities (v), MAC grid
        for (int j = 1; j < mResolution.y(); ++j)
        { // skip first and last column: those are determined by the boundary condition
            for (int i = 0; i < mResolution.x(); ++i)
            {
                // TODO: Compute the velocity
                float v_current = mVelocity_v->getVertexDataAt({ i, j }).x();
                float u_bottom = (mVelocity_u->getVertexDataAt({ i, j - 1 }).x() + 
                                  mVelocity_u->getVertexDataAt({ i + 1, j - 1 }).x()) / 2.0f;
                float u_top = (mVelocity_u->getVertexDataAt({ i, j }).x() + 
                                  mVelocity_u->getVertexDataAt({ i + 1, j }).x()) / 2.0f;
                
                float last_x_velocity = (u_bottom + u_top) / 2.0f;
                float last_y_velocity = v_current;

                // TODO: Find the last position of the particle (in grid coordinates) using an Euler step
                float last_x = i - stepSize * last_x_velocity / mSpacing.x();
                float last_y = j - stepSize * last_y_velocity / mSpacing.y();

                // Make sure the coordinates are inside the boundaries
                if (last_x < 0)
                    last_x = 0;
                if (last_y < 0)
                    last_y = 0;
                if (last_x > mResolution.x() - 1)
                    last_x = mResolution.x() - 1;
                if (last_y > mResolution.y() - 0)
                    last_y = mResolution.y() - 0;

                // Determine corners for bilinear interpolation
                int x_low  = (int)last_x;
                int y_low  = (int)last_y;
                int x_high = std::min(x_low + 1, mResolution.x() - 1);
                int y_high = std::min(y_low + 1, mResolution.y() - 0);

                // Compute the interpolation weights
                float x_weight = last_x - x_low;
                float y_weight = last_y - y_low;

                // TODO: Bilinear interpolation
                float v00 = mVelocity_v->getVertexDataAt({ x_low, y_low }).x();
                float v10 = mVelocity_v->getVertexDataAt({ x_high, y_low }).x();
                float v01 = mVelocity_v->getVertexDataAt({ x_low, y_high }).x();
                float v11 = mVelocity_v->getVertexDataAt({ x_high, y_high }).x();
                
                float interpolated_v = (1 - x_weight) * (1 - y_weight) * v00 +
                                      x_weight * (1 - y_weight) * v10 +
                                      (1 - x_weight) * y_weight * v01 +
                                      x_weight * y_weight * v11;
                
                mVelocity_v_tmp->setVertexDataAt({ i, j }, interpolated_v);
            }
        }

        // Copy the values in temp to the original buffers
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 1; i < mResolution.x(); ++i)
                mVelocity_u->setVertexDataAt({ i, j }, mVelocity_u_tmp->getVertexDataAt({ i, j }));
        for (int j = 1; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                mVelocity_v->setVertexDataAt({ i, j }, mVelocity_v_tmp->getVertexDataAt({ i, j }));
    }

    void FluidSolver::advectDensity()
    {
        // Densities live on the grid centers, the velocities on the MAC grid
        // Separate their computation to avoid confusion

        //     the staggered grid indices look like this (with half offsets added already)
        //     ------------------------------
        //     |                            |
        //     O u_i,j       X p_i,j        O u_i+1,j
        //     |                            |
        //     ------------------------------

        // Densities, grid centers
        for (int j = 0; j < mResolution.y(); ++j)
        {
            for (int i = 0; i < mResolution.x(); ++i)
            {
                // TODO: Compute the velocity
                float u_left = mVelocity_u->getVertexDataAt({ i, j }).x();
                float u_right = mVelocity_u->getVertexDataAt({ i + 1, j }).x();
                float v_bottom = mVelocity_v->getVertexDataAt({ i, j }).x();
                float v_top = mVelocity_v->getVertexDataAt({ i, j + 1 }).x();
                
                float last_x_velocity = (u_left + u_right) / 2.0f;
                float last_y_velocity = (v_bottom + v_top) / 2.0f;

                // TODO: Find the last position of the particle (in grid coordinates) using an Euler step
                float last_x = i - stepSize * last_x_velocity / mSpacing.x();
                float last_y = j - stepSize * last_y_velocity / mSpacing.y();

                // Make sure the coordinates are inside the boundaries
                const float offset = 0.0; // a trick to fight the dissipation through boundaries is to sample with a small offset
                if (last_x < offset)
                    last_x = offset;
                if (last_y < offset)
                    last_y = offset;
                if (last_x > mResolution.x() - 1 - offset)
                    last_x = mResolution.x() - 1 - offset;
                if (last_y > mResolution.y() - 1 - offset)
                    last_y = mResolution.y() - 1 - offset;

                // Determine corners for bilinear interpolation
                int x_low  = (int)last_x;
                int y_low  = (int)last_y;
                int x_high = std::min(x_low + 1, mResolution.x() - 1);
                int y_high = std::min(y_low + 1, mResolution.y() - 1);

                // Compute the interpolation weights
                float x_weight = last_x - x_low;
                float y_weight = last_y - y_low;

                // TODO: Bilinear interpolation
                float rho00 = mDensity->getVertexDataAt({ x_low, y_low }).x();
                float rho10 = mDensity->getVertexDataAt({ x_high, y_low }).x();
                float rho01 = mDensity->getVertexDataAt({ x_low, y_high }).x();
                float rho11 = mDensity->getVertexDataAt({ x_high, y_high }).x();
                
                float interpolated_rho = (1 - x_weight) * (1 - y_weight) * rho00 +
                                        x_weight * (1 - y_weight) * rho10 +
                                        (1 - x_weight) * y_weight * rho01 +
                                        x_weight * y_weight * rho11;
                
                mDensity_tmp->setVertexDataAt({ i, j }, interpolated_rho);
            }
        }

        // Copy the values in temp to the original buffers
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                mDensity->setVertexDataAt({ i, j }, mDensity_tmp->getVertexDataAt({ i, j }));
    }

    void FluidSolver::buildLaplace2d(Eigen::SparseMatrix<double>& A)
    {
        Eigen::Vector2i resolution = mPressure->getGrid()->getResolution();

        std::list<Eigen::Triplet<double>> coeff;
        for (int j = 0; j < resolution.y(); ++j)
            for (int i = 0; i < resolution.x(); ++i)
            {

                if (i > 0)
                    coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i - 1, 1));
                if (i < resolution.x() - 1)
                    coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i + 1, 1));
                if (j > 0)
                    coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, (j - 1) * resolution.x() + i, 1));
                if (j < resolution.y() - 1)
                    coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, (j + 1) * resolution.x() + i, 1));

                switch (boundary)
                {
                case EBoundary::Open: // dirichlet condition: smooth continuation
                    // example for 4x4 matrix
                    // -4  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0
                    //  1 -4  1  0  0  1  0  0  0  0  0  0  0  0  0  0
                    // 	0  1 -4  1  0  0  1  0  0  0  0  0  0  0  0  0
                    // 	0  0  1 -4  0  0  0  1  0  0  0  0  0  0  0  0
                    //  1  0  0  0 -4  1  0  0  1  0  0  0  0  0  0  0
                    // 	0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0  0
                    // 	0  0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0
                    // 	0  0  0  1  0  0  1 -4  0  0  0  1  0  0  0  0
                    // 	0  0  0  0  1  0  0  0 -4  1  0  0  1  0  0  0
                    // 	0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0  0
                    // 	0  0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0
                    // 	0  0  0  0  0  0  0  1  0  0  1 -4  0  0  0  1
                    // 	0  0  0  0  0  0  0  0  1  0  0  0 -4  1  0  0
                    // 	0  0  0  0  0  0  0  0  0  1  0  0  1 -4  1  0
                    // 	0  0  0  0  0  0  0  0  0  0  1  0  0  1 -4  1
                    // 	0  0  0  0  0  0  0  0  0  0  0  1  0  0  1 -4

                    coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i, -4));

                    break;

                case EBoundary::Closed: // Neumann condition: forward/backward difference on boundary   (same pattern of -1's, but the main diagonal contains the number of -1's per row)
                    // -2  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0
                    //  1 -3  1  0  0  1  0  0  0  0  0  0  0  0  0  0
                    //	0  1 -3  1  0  0  1  0  0  0  0  0  0  0  0  0
                    //	0  0  1 -2  0  0  0  1  0  0  0  0  0  0  0  0
                    //  1  0  0  0 -3  1  0  0  1  0  0  0  0  0  0  0
                    //	0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0  0
                    //	0  0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0
                    //	0  0  0  1  0  0  1 -3  0  0  0  1  0  0  0  0
                    //	0  0  0  0  1  0  0  0 -3  1  0  0  1  0  0  0
                    //	0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0  0
                    //	0  0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0
                    //	0  0  0  0  0  0  0  1  0  0  1 -3  0  0  0  1
                    //	0  0  0  0  0  0  0  0  1  0  0  0 -2  1  0  0
                    //	0  0  0  0  0  0  0  0  0  1  0  0  1 -3  1  0
                    //	0  0  0  0  0  0  0  0  0  0  1  0  0  1 -3  1
                    //	0  0  0  0  0  0  0  0  0  0  0  1  0  0  1 -2

                    int k= 0; // number of -1's per row
                    if (i > 0) k++;
                    if (i < resolution.x() - 1) k++;
                    if (j > 0) k++;
                    if (j < resolution.y() - 1) k++;
                    coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i, -k));

                    break;
                }
            }

        A.setFromTriplets(coeff.begin(), coeff.end());
    }
}
