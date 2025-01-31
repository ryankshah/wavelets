defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform
  """

  alias Wavelets.{DWT, Filter, IDWT}

  @doc """
  Performs n-dimensional forward discrete wavelet transform.
  The input signal should be a nested list structure matching the dimensions.
  Returns a tuple with approximation and a list of detail coefficients.
  """
  @spec forward(list(), Filter.t(), non_neg_integer(), Wavelets.precision()) ::
          {list(), list()}
  def forward(signal, filter, dims, _precision \\ :double) do
    case dims do
      1 -> DWT.forward_1d(signal, filter)
      2 -> DWT.forward_2d(signal, filter)
      _ -> do_forward_nd(signal, filter, dims)
    end
  end

  @doc """
  Performs n-dimensional inverse discrete wavelet transform
  """
  @spec inverse(
          list(),
          list(),
          Filter.t(),
          non_neg_integer(),
          Wavelets.precision()
        ) :: list()
  def inverse(approximation, details, filter, dims, _precision \\ :double) do
    case dims do
      1 -> IDWT.inverse_1d(approximation, details, filter)
      2 -> IDWT.inverse_2d(approximation, details, filter)
      _ -> do_inverse_nd(approximation, details, filter, dims)
    end
  end

  # Private helper functions
  defp do_forward_nd(signal, filter, dims) do
    coeffs = apply_forward_transforms(signal, filter, dims)
    extract_coeffs(coeffs, dims)
  end

  defp apply_forward_transforms(signal, filter, dims) do
    Enum.reduce((dims - 1)..0, signal, fn dim, acc ->
      transform_dimension(acc, dim, filter)
    end)
  end

  defp transform_dimension(signal_nd, 0, filter) do
    # Base case: apply 1D transform to the innermost lists
    Enum.map(signal_nd, &DWT.forward_1d(&1, filter))
  end

  defp transform_dimension(signal_nd, dim, filter) do
    # Recursive case: traverse to inner dimensions
    Enum.map(signal_nd, &transform_dimension(&1, dim - 1, filter))
  end

  defp do_inverse_nd(approximation, details, filter, dims) do
    coeffs = merge_coeffs(approximation, details, dims)
    apply_inverse_transforms(coeffs, filter, dims)
  end

  defp apply_inverse_transforms(coeffs, filter, dims) do
    Enum.reduce(0..(dims - 1), coeffs, fn dim, acc ->
      inverse_dimension(acc, dim, filter)
    end)
  end

  defp inverse_dimension(signal_nd, 0, filter) do
    # Base case: apply 1D inverse transform
    Enum.map(signal_nd, fn {a, d} -> IDWT.inverse_1d(a, d, filter) end)
  end

  defp inverse_dimension(signal_nd, dim, filter) do
    # Recursive case: traverse to inner dimensions
    Enum.map(signal_nd, &inverse_dimension(&1, dim - 1, filter))
  end

  defp extract_coeffs(coeffs, 1), do: Enum.unzip(coeffs)

  defp extract_coeffs(coeffs, dims) do
    {approx, details} =
      coeffs
      |> Enum.map(&extract_coeffs(&1, dims - 1))
      |> Enum.unzip()

    {approx, List.flatten(details)}
  end

  defp merge_coeffs(approximation, details, 1) do
    Enum.zip(approximation, details)
  end

  defp merge_coeffs(approximation, details, dims) do
    detail_chunks = Enum.chunk_every(details, div(length(details), dims))

    Enum.zip(approximation, detail_chunks)
    |> Enum.map(fn {a, d} -> merge_coeffs(a, d, dims - 1) end)
  end
end
